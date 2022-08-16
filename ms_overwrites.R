library(measurements)
ms_conversions <- function(d,
                           convert_units_from = 'mg/l',
                           convert_units_to,
                           convert_molecules,
                           macrosheds_root){

    if(missing(macrosheds_root)){
        stop('Please provide macrosheds_root, information needed to convert variables is stored here')
    }

    ms_vars_path <- paste0(macrosheds_root, '/ms_vars.feather')

    if(! file.exists(ms_vars_path)){
        ms_vars <- readr::read_csv('https://figshare.com/articles/dataset/variable_metadata/19358585/files/35134504',
                                   col_types = readr::cols())

        feather::write_feather(ms_vars, ms_vars_path)
    } else{
        ms_vars <- feather::read_feather(ms_vars_path)
    }

    #checks
    cm <- ! missing(convert_molecules)
    cuF <- ! missing(convert_units_from) && ! is.null(convert_units_from)
    cuT <- ! missing(convert_units_to) && ! is.null(convert_units_to)

    if(sum(cuF, cuT) == 1){
        stop('convert_units_from and convert_units_to must be supplied together')
    }
    if(length(convert_units_from) != length(convert_units_to)){
        stop('convert_units_from and convert_units_to must have the same length')
    }

    vars <- ms_drop_var_prefix(d$var)

    if(any(!vars %in% ms_vars$variable_code)){
        not_a_ms_var <- unique(vars[!vars %in% ms_vars$variable_code])
        stop(paste0(paste(not_a_ms_var, collapse = ', '),
                    ' is not a MacroSheds variable. only MacroSheds variables can be converted'))
    }

    if(any(duplicated(names(convert_units_from)))){
        stop('duplicated names in convert_units_from')
    }
    if(any(duplicated(names(convert_units_to)))){
        stop('duplicated names in convert_units_to')
    }

    vars_convertable <- ms_vars %>%
        filter(variable_code %in% !!vars) %>%
        pull(unit) %>%
        tolower()

    if(length(convert_units_from) == 1 && length(convert_units_to) == 1){
        if(! all(vars_convertable == 'mg/l')){
            print(all(vars_convertable))
            warning('unable to convert non-concentration variables')
        }
    } else{
        if(! all(vars %in% names(convert_units_from)) || ! all(vars %in% names(convert_units_to))){
            stop('when specifying individual variable conversions, all variables in d must be accounted for')
        }
            cu_shared_names <- base::intersect(names(convert_units_from),
                                               names(convert_units_to))

            if(length(cu_shared_names) != length(convert_units_to)){
                stop('names of convert_units_from and convert_units_to must match')
            }
    }

    convert_units_from <- tolower(convert_units_from)
    convert_units_to <- tolower(convert_units_to)

    whole_molecule <- c('NO3', 'SO4', 'PO4', 'SiO2', 'SiO3', 'NH4', 'NH3',
                        'NO3_NO2')
    element_molecule <- c('NO3_N', 'SO4_S', 'PO4_P', 'SiO2_S', 'SiO3_S', 'NH4_N',
                          'NH3_N', 'NO3_NO2_N')

    if(cm){
        whole_to_element <- grep(paste0(paste0('^', convert_molecules, '$'), collapse = '|'),
                                 whole_molecule)
        element_to_whole <- grep(paste0(paste0('^', convert_molecules, '$'), collapse = '|'),
                                 element_molecule)

        if(length(element_to_whole) == 0 && length(whole_to_element) == 0){
            stop(paste0('convert_molecules must be one of: ', paste(whole_molecule, collapse = ' '),
                        ' or: ', paste(element_molecule, collapse = ' ')))
        }
    } else{
        convert_molecules <- NULL
    }

    molecular_conversion_map <- list(
        NH4 = 'N',
        NO3 = 'N',
        NH3 = 'N',
        SiO2 = 'Si',
        SiO3 = 'Si',
        SO4 = 'S',
        PO4 = 'P',
        NO3_NO2 = 'N')

    # handle molecular conversions, like NO3 -> NO3_N
    if(cm && length(whole_to_element) > 0){
        convert_molecules_element <-  whole_molecule[whole_to_element]
        for(v in 1:length(convert_molecules_element)){

            molecule_real <- ms_vars %>%
                filter(variable_code == !!convert_molecules_element[v]) %>%
                pull(molecule)

            if(is.na(molecule_real)) {
                molecule_real <- convert_molecules_element[v]
            }

            d$val[vars == convert_molecules_element[v]] <-
                convert_molecule(x = d$val[vars == convert_molecules_element[v]],
                                 from = molecule_real,
                                 to = unname(molecular_conversion_map[v]))

            check_double <- stringr::str_split_fixed(unname(molecular_conversion_map[v]), '', n = Inf)[1,]

            if(length(check_double) > 1 && length(unique(check_double)) == 1) {
                molecular_conversion_map[v] <- unique(check_double)
            }

            new_name <- paste0(d$var[vars == convert_molecules_element[v]], '_', unname(molecular_conversion_map[v]))

            d$var[vars == convert_molecules_element[v]] <- new_name
        }
    }

    # handle molecular conversions, like NO3_N -> NO3
    if(cm && length(element_to_whole) > 0){
        convert_molecules_element <-  element_molecule[element_to_whole]
        for(v in 1:length(convert_molecules_element)){

            molecule_real <- ms_vars %>%
                filter(variable_code == !!convert_molecules_element[v]) %>%
                pull(molecule)

            if(is.na(molecule_real)) {
                molecule_real <- convert_molecules_element[v]
            }

            d$val[vars == convert_molecules_element[v]] <-
                convert_molecule(x = d$val[vars == convert_molecules_element[v]],
                                 from = molecule_real,
                                 to = whole_molecule[element_to_whole[v]])

            # check_double <- stringr::str_split_fixed(unname(molecular_conversion_map[v]), '', n = Inf)[1,]
            #
            # if(length(check_double) > 1 && length(unique(check_double)) == 1) {
            #     molecular_conversion_map[v] <- unique(check_double)
            # }
            old_var <- unique(d$var[vars == convert_molecules_element[v]])
            new_name <- substr(d$var[vars == convert_molecules_element[v]], 0, nchar(old_var)-2)

            d$var[vars == convert_molecules_element[v]] <- new_name
        }
    }

    # Turn a single input into a named vector with all variables in dataframe
    if(length(convert_units_from) == 1){
        all_vars <- unique(vars)
        convert_units_from <- rep(convert_units_from, length(all_vars))
        names(convert_units_from) <- all_vars
        convert_units_to <- rep(convert_units_to, length(all_vars))
        names(convert_units_to) <- all_vars
    }

    # Converts input to grams if the final unit contains grams
    for(i in 1:length(convert_units_from)){

        unitfrom <- convert_units_from[i]
        unitto <- convert_units_to[i]
        v <- names(unitfrom)

        g_conver <- FALSE
        if(grepl('mol|eq', unitfrom) && grepl('g', unitto) || v %in% convert_molecules){

            molecule_real <- ms_vars %>%
                filter(variable_code == !!v) %>%
                pull(molecule)

            if(! is.na(molecule_real)){
                formula <- molecule_real
            } else {
                formula <- v
            }

            d$val[vars == v] <- convert_to_gl(x = d$val[vars == v],
                                              input_unit = unitfrom,
                                              formula = formula,
                                              ms_vars = ms_vars)

            g_conver <- TRUE
        }

        #convert prefix
        print('UNITS OCNVERT MAGIC TIME')
        print(unitfrom)
        print(unitfrom[[1]])
        print(unitto)
        print(unitto[[1]])


        d$val[vars == v] <- convert_unit(x = d$val[vars == v],
                                         input_unit = unitfrom,
                                         output_unit = unitto)

        d$val[vars == v] <- convert_unit(x = d$val[vars == v],
                                         input_unit = unitfrom,
                                         output_unit = unitto)

        #Convert to mol or eq if that is the output unit
        if(grepl('mol|eq', unitto)) {

            d$val[vars == v] <- convert_from_gl(x = d$val[vars == v],
                                                input_unit = unitfrom,
                                                output_unit = unitto,
                                                molecule = v,
                                                g_conver = g_conver,
                                                ms_vars = ms_vars)
        }
    }

    return(d)
}
convert_unit <- function(x, input_unit, output_unit){

    units <- tibble(prefix = c('n', "u", "m", "c", "d", "h", "k", "M"),
                    convert_factor = c(0.000000001, 0.000001, 0.001, 0.01, 0.1, 100,
                                       1000, 1000000))

    old_fraction <- as.vector(stringr::str_split_fixed(input_unit, "/", n = Inf))
    old_top <- as.vector(stringr::str_split_fixed(old_fraction[1], "", n = Inf))

    if(length(old_fraction) == 2) {
        old_bottom <- as.vector(stringr::str_split_fixed(old_fraction[2], "", n = Inf))
    }

    new_fraction <- as.vector(stringr::str_split_fixed(output_unit, "/", n = Inf))
    new_top <- as.vector(stringr::str_split_fixed(new_fraction[1], "", n = Inf))

    if(length(new_fraction == 2)) {
        new_bottom <- as.vector(stringr::str_split_fixed(new_fraction[2], "", n = Inf))
    }

    old_top_unit <- tolower(stringr::str_split_fixed(old_top, "", 2)[1])

    if(old_top_unit %in% c('g', 'e', 'q', 'l') || old_fraction[1] == 'mol') {
        old_top_conver <- 1
    } else {
        old_top_conver <- as.numeric(filter(units, prefix == old_top_unit)[,2])
    }

    old_bottom_unit <- tolower(stringr::str_split_fixed(old_bottom, "", 2)[1])

    if(old_bottom_unit %in% c('g', 'e', 'q', 'l') || old_fraction[2] == 'mol') {
        old_bottom_conver <- 1
    } else {
        old_bottom_conver <- as.numeric(filter(units, prefix == old_bottom_unit)[,2])
    }

    new_top_unit <- tolower(stringr::str_split_fixed(new_top, "", 2)[1])

    if(new_top_unit %in% c('g', 'e', 'q', 'l') || new_fraction[1] == 'mol') {
        new_top_conver <- 1
    } else {
        new_top_conver <- as.numeric(filter(units, prefix == new_top_unit)[,2])
    }

    new_bottom_unit <- tolower(stringr::str_split_fixed(new_bottom, "", 2)[1])

    if(new_bottom_unit %in% c('g', 'e', 'q', 'l') || new_fraction[2] == 'mol') {
        new_bottom_conver <- 1
    } else {
        new_bottom_conver <- as.numeric(filter(units, prefix == new_bottom_unit)[,2])
    }

    new_val <- x*old_top_conver
    new_val <- new_val/new_top_conver

    new_val <- new_val/old_bottom_conver
    new_val <- new_val*new_bottom_conver

    return(new_val)
}

# End unit converstion
sd_or_0 <- function(x, na.rm = FALSE) {

    #Only used to bypass the tyranny of the errors package not letting
    #me take the mean of an errors object of length 1 without setting the
    #uncertainty to 0

    x <- if(is.vector(x) || is.factor(x)) x else as.double(x)

    if(length(x) == 1) return(0)

    x <- sqrt(var(x, na.rm = na.rm))
}

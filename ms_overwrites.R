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

# 'errors' package handlers
sd_or_0 <- function(x, na.rm = FALSE) {

    #Only used to bypass the tyranny of the errors package not letting
    #me take the mean of an errors object of length 1 without setting the
    #uncertainty to 0

    x <- if(is.vector(x) || is.factor(x)) x else as.double(x)

    if(length(x) == 1) return(0)

    x <- sqrt(var(x, na.rm = na.rm))
}

mean_or_x <- function(x, na.rm = FALSE) {
    # also used to bypass the tyranny of the errors package not letting
    # someone take the mean of an errors object of length 1. this func returns the
    # original value if the group is length one, and the mean otherwise

    if(length(x) == 1) return(x)

    x <- mean(var(x, na.rm = na.rm))
    print('multiple values meaned')
    return(x)
}


Mode <- function(x, na.rm = TRUE){

    if(na.rm){
        x <- na.omit(x)
    }

    ux <- unique(x)
    mode_out <- ux[which.max(tabulate(match(x, ux)))]
    return(mode_out)

}


approxjoin_datetime <- function(x,
                                y,
                                rollmax = '7:30',
                                keep_datetimes_from = 'x',
                                indices_only = FALSE){
    #direction = 'forward'){

    #x and y: macrosheds standard tibbles with only one site_code,
    #   which must be the same in x and y. Nonstandard tibbles may also work,
    #   so long as they have datetime columns, but the only case where we need
    #   this for other tibbles is inside precip_pchem_pflux_idw, in which case
    #   indices_only == TRUE, so it's not really set up for general-purpose joining
    #rollmax: the maximum snap time for matching elements of x and y.
    #   either '7:30' for continuous data or '12:00:00' for grab data
    #direction [REMOVED]: either 'forward', meaning elements of x will be rolled forward
    #   in time to match the next y, or 'backward', meaning elements of
    #   x will be rolled back in time to reach the previous y
    #keep_datetimes_from: string. either 'x' or 'y'. the datetime column from
    #   the corresponding tibble will be kept, and the other will be dropped
    #indices_only: logical. if TRUE, a join is not performed. rather,
    #   the matching indices from each tibble are returned as a named list of vectors..

    #good datasets for testing this function:
    # x <- tribble(
    #     ~datetime, ~site_code, ~var, ~val, ~ms_status, ~ms_interp,
    #     '1968-10-09 04:42:00', 'GSWS10', 'GN_alk', set_errors(27.75, 1), 0, 0,
    #     '1968-10-09 04:44:00', 'GSWS10', 'GN_alk', set_errors(21.29, 1), 0, 0,
    #     '1968-10-09 04:47:00', 'GSWS10', 'GN_alk', set_errors(21.29, 1), 0, 0,
    #     '1968-10-09 04:59:59', 'GSWS10', 'GN_alk', set_errors(16.04, 1), 0, 0,
    #     '1968-10-09 05:15:01', 'GSWS10', 'GN_alk', set_errors(17.21, 1), 1, 0,
    #     '1968-10-09 05:30:59', 'GSWS10', 'GN_alk', set_errors(16.50, 1), 0, 0) %>%
    # mutate(datetime = as.POSIXct(datetime, tz = 'UTC'))
    # y <- tribble(
    #     ~datetime, ~site_code, ~var, ~val, ~ms_status, ~ms_interp,
    #     '1968-10-09 04:00:00', 'GSWS10', 'GN_alk', set_errors(1.009, 1), 1, 0,
    #     '1968-10-09 04:15:00', 'GSWS10', 'GN_alk', set_errors(2.009, 1), 1, 1,
    #     '1968-10-09 04:30:00', 'GSWS10', 'GN_alk', set_errors(3.009, 1), 1, 1,
    #     '1968-10-09 04:45:00', 'GSWS10', 'GN_alk', set_errors(4.009, 1), 1, 1,
    #     '1968-10-09 05:00:00', 'GSWS10', 'GN_alk', set_errors(5.009, 1), 1, 1,
    #     '1968-10-09 05:15:00', 'GSWS10', 'GN_alk', set_errors(6.009, 1), 1, 1) %>%
    #     mutate(datetime = as.POSIXct(datetime, tz = 'UTC'))

    #tests
    if('site_code' %in% colnames(x) && length(unique(x$site_code)) > 1){
        stop('Only one site_code allowed in x at the moment')
    }
    if('var' %in% colnames(x) && length(unique(drop_var_prefix(x$var))) > 1){
        stop('Only one var allowed in x at the moment (not including prefix)')
    }
    if('site_code' %in% colnames(y) && length(unique(y$site_code)) > 1){
        stop('Only one site_code allowed in y at the moment')
    }
    if('var' %in% colnames(y) && length(unique(drop_var_prefix(y$var))) > 1){
        stop('Only one var allowed in y at the moment (not including prefix)')
    }
    if('site_code' %in% colnames(x) &&
       'site_code' %in% colnames(y) &&
       x$site_code[1] != y$site_code[1]) stop('x and y site_code must be the same')
    if(! rollmax %in% c('7:30', '12:00:00')) stop('rollmax must be "7:30" or "12:00:00"')
    # if(! direction %in% c('forward', 'backward')) stop('direction must be "forward" or "backward"')
    if(! keep_datetimes_from %in% c('x', 'y')) stop('keep_datetimes_from must be "x" or "y"')
    if(! 'datetime' %in% colnames(x) || ! 'datetime' %in% colnames(y)){
        stop('both x and y must have "datetime" columns containing POSIXct values')
    }
    if(! is.logical(indices_only)) stop('indices_only must be a logical')

    #deal with the case of x or y being a specialized "flow" tibble
    # x_is_flowtibble <- y_is_flowtibble <- FALSE
    # if('flow' %in% colnames(x)) x_is_flowtibble <- TRUE
    # if('flow' %in% colnames(y)) y_is_flowtibble <- TRUE
    # if(x_is_flowtibble && ! y_is_flowtibble){
    #     varname <- y$var[1]
    #     y$var = NULL
    # } else if(y_is_flowtibble && ! x_is_flowtibble){
    #     varname <- x$var[1]
    #     x$var = NULL
    # } else if(! x_is_flowtibble && ! y_is_flowtibble){
    #     varname <- x$var[1]
    #     x$var = NULL
    #     y$var = NULL
    # } else {
    #     stop('x and y are both "flow" tibbles. There should be no need for this')
    # }
    # if(x_is_flowtibble) x <- rename(x, val = flow)
    # if(y_is_flowtibble) y <- rename(y, val = flow)

    #data.table doesn't work with the errors package, so error needs
    #to be separated into its own column. also give same-name columns suffixes

    if('val' %in% colnames(x)){ #crude catch for nonstandard ms tibbles (fine for now)
        x <- x %>%
            mutate(err = errors::errors(val),
                   val = errors::drop_errors(val)) %>%
            rename_with(.fn = ~paste0(., '_x'),
                        .cols = everything()) %>%
            # .cols = any_of(c('site_code', 'var', 'val',
            #                  'ms_status', 'ms_interp'))) %>%
            data.table::as.data.table()

        y <- y %>%
            mutate(err = errors::errors(val),
                   val = errors::drop_errors(val)) %>%
            rename_with(.fn = ~paste0(., '_y'),
                        .cols = everything()) %>%
            data.table::as.data.table()
    } else {
        x <- dplyr::rename(x, datetime_x = datetime) %>% data.table::as.data.table()
        y <- dplyr::rename(y, datetime_y = datetime) %>% data.table::as.data.table()
    }

    #alternative implementation of the "on" argument in data.table joins...
    #probably more flexible, so leaving it here in case we need to do something crazy
    # data.table::setkeyv(x, 'datetime')
    # data.table::setkeyv(y, 'datetime')

    #convert the desired maximum roll distance from string to integer seconds
    rollmax <- ifelse(test = rollmax == '7:30',
                      yes = 7 * 60 + 30,
                      no = 12 * 60 * 60)

    #leaving this here in case the nearest neighbor join implemented below is too
    #slow. then we can fall back to a basic rolling join with a maximum distance
    # rollmax <- ifelse(test = direction == 'forward',
    #                   yes = -rollmax,
    #                   no = rollmax)
    #rollends will move the first/last value of x in the opposite `direction` if necessary
    # joined <- y[x, on = 'datetime', roll = rollmax, rollends = c(TRUE, TRUE)]

    #create columns in x that represent the snapping window around each datetime
    x[, `:=` (datetime_min = datetime_x - rollmax,
              datetime_max = datetime_x + rollmax)]
    y[, `:=` (datetime_y_orig = datetime_y)] #datetime col will be dropped from y

    # if(indices_only){
    #     y_indices <- y[x,
    #                    on = .(datetime_y <= datetime_max,
    #                           datetime_y >= datetime_min),
    #                    which = TRUE]
    #     return(y_indices)
    # }

    #join x rows to y if y's datetime falls within the x range
    joined <- y[x, on = .(datetime_y <= datetime_max,
                          datetime_y >= datetime_min)]
    joined <- na.omit(joined, cols = 'datetime_y_orig') #drop rows without matches

    #for any datetimes in x or y that were matched more than once, keep only
    #the nearest match
    joined[, `:=` (datetime_match_diff = abs(datetime_x - datetime_y_orig))]
    joined <- joined[, .SD[which.min(datetime_match_diff)], by = datetime_x]
    joined <- joined[, .SD[which.min(datetime_match_diff)], by = datetime_y_orig]

    if(indices_only){
        y_indices <- which(y$datetime_y %in% joined$datetime_y_orig)
        x_indices <- which(x$datetime_x %in% joined$datetime_x)
        return(list(x = x_indices, y = y_indices))
    }

    #drop and rename columns (data.table makes weird name modifications)
    if(keep_datetimes_from == 'x'){
        joined[, c('datetime_y', 'datetime_y.1', 'datetime_y_orig', 'datetime_match_diff') := NULL]
        data.table::setnames(joined, 'datetime_x', 'datetime')
    } else {
        joined[, c('datetime_x', 'datetime_y.1', 'datetime_y', 'datetime_match_diff') := NULL]
        data.table::setnames(joined, 'datetime_y_orig', 'datetime')
    }

    #restore error objects, var column, original column names (with suffixes).
    #original column order
    joined <- as_tibble(joined) %>%
        mutate(val_x = errors::set_errors(val_x, err_x),
               val_y = errors::set_errors(val_y, err_y)) %>%
        select(-err_x, -err_y)
    # mutate(var = !!varname)

    # if(x_is_flowtibble) joined <- rename(joined,
    #                                      flow = val_x,
    #                                      ms_status_flow = ms_status_x,
    #                                      ms_interp_flow = ms_interp_x)
    # if(y_is_flowtibble) joined <- rename(joined,
    #                                      flow = val_y,
    #                                      ms_status_flow = ms_status_y,
    #                                      ms_interp_flow = ms_interp_y)

    # if(! sum(grepl('^val_[xy]$', colnames(joined))) > 1){
    #     joined <- rename(joined, val = matches('^val_[xy]$'))
    # }

    joined <- select(joined,
                     datetime,
                     # matches('^val_?[xy]?$'),
                     # any_of('flow'),
                     starts_with('site_code'),
                     any_of(c(starts_with('var_'), matches('^var$'))),
                     any_of(c(starts_with('val_'), matches('^val$'))),
                     starts_with('ms_status_'),
                     starts_with('ms_interp_'))

    return(joined)
}

ms_calc_flux <- function(chemistry, q, q_type, site_info = NULL, verbose = TRUE,
                         method = 'simple', aggregation = 'simple') {

    #### Checks
    if(! all(c('site_code', 'val', 'var', 'datetime', 'ms_interp', 'ms_status') %in% names(chemistry))){
        stop('The argument to chemistry must contain precipitation chemistry or stream chemistry data in MacroSheds format (column names of site_code, val, var, datetime, ms_interp, ms_status at minimum).')
    }
    if(! all(c('site_code', 'val', 'var', 'datetime', 'ms_interp', 'ms_status') %in% names(q))){
        stop('The argument to q must contain precipitation or stream discharge data in MacroSheds format (column names of site_code, val, var, datetime, ms_interp, ms_status at minimum).')
    }
    if(! grepl('(precipitation|discharge)', q_type)){
        stop('q_type must be "discharge" or "precipitation"')
    }
    if(! 'POSIXct' %in% class(q$datetime)){
        q$datetime <- as.POSIXct(q$datetime)
    }
    if(! 'POSIXct' %in% class(chemistry$datetime)){
        chemistry$datetime <- as.POSIXct(chemistry$datetime)
    }

    # check that method, if non-null, is in accepted list
    rsfme_accepted <- c('average', 'pw', 'composite', 'wrtds', 'beale', 'simple')
    if(!method %in% rsfme_accepted) {
      stop(glue('method supplied is not in accepted list, must be one of the following:\n {list}',
                list = rsfme_accepted))
    } else {
      writeLines(glue('calculating flux using method: {method}', method = method))
    }

    # make sure agg option is annual or monthly if calculating any non-null method
    # and otherwise timestep is data-res and using simple QC
    rsfme_aggs <- c('annual', 'monthly', 'simple')
    if(!aggregation %in% rsfme_aggs) {
      stop(glue('time aggregation is not in accepted list, must be one of the following:\n {list}',
                list = rsfme_aggs))
    } else if(aggregation == 'simple') {
      writeLines(glue('calculating flux at highest possible resolution timestep of data supplied, using simple Q*C methods', aggregation = aggregation))
    } else {
      writeLines(glue('calculating flux over: {aggregation}', aggregation = aggregation))
    }

    if(q_type == 'discharge' && is.null(site_info)) {
        site_info <- try(ms_download_site_data())

        if(inherits(site_info, 'try-error')){
            stop("When q_type == 'discharge', you must either have site_info defined as the MacroSheds \n
                 site_data table or you must have an internet connection to download the table with ms_download_site_data()")
        }
    } else {
        site_info$ws_area_ha <- errors::set_errors(site_info$ws_area_ha, 0)
    }

    # Check both files have the same sites
    sites_chem <- unique(chemistry$site_code)
    sites_q <- unique(q$site_code)

    if(! all(sites_chem %in% sites_q)){
        stop('Both chemistry and q must have the same sites')
    }

    sites <- sites_chem

    # Check the intervals are the same in both chemistry and q
    q_interval <- Mode(diff(as.numeric(q$datetime)))

    interval <- case_when(q_interval == 86400 ~ 'daily',
                          q_interval == 3600 ~ 'hourly',
                          q_interval == 1800 ~ '30 minute',
                          q_interval == 960 ~ '15 minute',
                          q_interval == 600 ~ '10 minute',
                          q_interval == 300 ~ '5 minute',
                          q_interval == 60 ~ '1 minute')

    # q_interval <- errors::as.errors(q_interval)


    flow_is_highres <- Mode(diff(as.numeric(q$datetime))) <= 15 * 60
    if(is.na(flow_is_highres)) { flow_is_highres <- FALSE }

    if(is.na(interval)) {
        stop(paste0('interval of samples must be one',
                    ' of: daily, hourly, 30 minute, 15 minute, 10 minute, 5 minute, or 1 minute.',
                    ' See macrosheds::ms_synchronize_timestep() to standardize your intervals.'))
    } else if(verbose) {
        print(paste0('q dataset has a ', interval, ' interval'))
    }

    # add errors if they don't exist
    if('val_err' %in% names(chemistry)){
        errors::errors(chemistry$val) <- chemistry$val_err

        chemistry <- chemistry %>%
            select(-val_err)

    } else if(all(errors::errors(chemistry$val) == 0)){
        errors::errors(chemistry$val) <- 0
    }

    if('val_err' %in% names(q)){
        errors::errors(q$val) <- q$val_err

        q <- q %>%
            select(-val_err)

    } else if(all(errors::errors(q$val) == 0)){
        errors::errors(q$val) <- 0
    }

    # calc flux
    all_sites_flux <- tibble()

    for(s in 1:length(sites)) {

        site <- sites[s]

        site_chem <- chemistry %>%
            filter(site_code == !!site)

        site_q <- q %>%
            filter(site_code == !!site)

        daterange <- range(site_chem$datetime)

        site_q <- site_q %>%
            filter(
                site_code == !!site,
                datetime >= !!daterange[1],
                datetime <= !!daterange[2])

        if(nrow(site_q) == 0) { return(NULL) }

        chem_split <- site_chem %>%
            group_by(var) %>%
            arrange(datetime) %>%
            dplyr::group_split() %>%
            as.list()

        # Loop though all variables
        for(i in 1:length(chem_split)) {
            # df of just one solute chem at one site, over all time
            chem_chunk <- chem_split[[i]]

            # target solute
            target_solute <- unique(chem_chunk %>% pull(var))

            writeLines(glue('________\n\nformula: {method}\nsolute: {solute}\n________', method = method, solute = target_solute ))

            # 'good year' checks for RSFME calcs
            if(method != 'simple') {

              # df to populate with annual flux values by method
              out_frame <- tibble(wy = as.character(),
                    site_code = as.character(),
                    val = as.numeric(),
                    var = as.character(),
                    method = as.character())
                    ## ms_reccomended = as.integer(),
                    ## ms_interp_ratio = as.numeric(),
                    ## ms_status_ratio = as.numeric(),
                    ## ms_missing_ratio = as.numeric())


              # find acceptable years
              q_check <- raw_data_q %>%
                  mutate(date = date(datetime)) %>%
                  # NOTE: should we filter out NAs?
                  filter(ms_interp == 0, !is.na(val)) %>%
                  distinct(., date, .keep_all = TRUE) %>%
                  mutate(water_year = water_year(datetime, origin = "usgs")) %>%
                  group_by(water_year) %>%
                  summarise(n = n()) %>%
                  filter(n >= 311)

              conc_check <- raw_data_con %>%
                  mutate(date = date(datetime)) %>%
                  # NOTE: should we filter out NAs?
                  filter(!is.na(val)) %>%
                  distinct(., date, .keep_all = TRUE) %>%
                  mutate(water_year = water_year(date, origin = "usgs"),
                         quart = quarter(date)) %>%
                  group_by(water_year) %>%
                  summarise(count = n_distinct(quart),
                            n = n()) %>%
                  filter(n >= 4,
                         count > 3)


              q_good_years <- q_check$water_year
              conc_good_years <- conc_check$water_year

              # 'good years' where Q and Chem data both meet min requirements
              good_years <- q_good_years[q_good_years %in% conc_good_years]
              n_yrs <- length(good_years)

              # NOTE: adding handling if concentration data fails conc check
              if(nrow(conc_check) < 1) {
                writeLines(glue("{site} concentration data insufficient sample size and frequency to warrant flux estimation",
                                "\n   no water years in {site} dataset with minimum standards met", site = site_code))
                next
              } else if(nrow(q_check) < 1) {
                writeLines(glue("{site} discharge data insufficient sample size and frequency to warrant flux estimation",
                                "\n   no water years in {site} dataset with minimum standards met", site = site_code))
                next
              } else if(length(good_years) == 0) {
                writeLines(glue("no water years where q data and concentration data both meet minimum standards",
                      "skipping site: {site}", site = site_code))
                next
              }

              #join data and cut to good years
              daily_data_con <- raw_data_con %>%
                  mutate(date = date(datetime)) %>%
                  group_by(date) %>%
                  summarize(val = mean_or_x(val)) %>%
                  mutate(site_code = !!site_code, var = 'con') %>%
                  select(site_code, datetime = date, var, val)

              daily_data_q <- raw_data_q %>%
                  mutate(date = date(datetime)) %>%
                  group_by(date) %>%
                  summarize(val = mean_or_x(val)) %>%
                  mutate(site_code = !!site_code, var = 'q_lps') %>%
                  select(site_code, datetime = date, var, val)

              q_df <- daily_data_q %>%
                pivot_wider(names_from = var,
                            values_from = val)

              raw_data_full <- rbind(daily_data_con, daily_data_q) %>%
                  pivot_wider(names_from = var, values_from = val, id_cols = c(site_code, datetime)) %>%
                 mutate(wy = water_year(datetime, origin = 'usgs')) %>%
                 filter(wy %in% good_years)

               con_full <- raw_data_full %>%
                   mutate(wy = as.numeric(as.character(wy))) %>%
                     select(site_code, datetime, con, wy) %>%
                     ## filter(wy < 1975) %>%
                 na.omit()

              if(tolower(method) == 'wrtds') {

                #### calculate WRTDS ######
                tryCatch(
                  expr = {
                    flux_annual_wrtds <- calculate_wrtds(
                      chem_df = con_full,
                      q_df = q_df,
                      ws_size = area,
                      lat = lat,
                      long = long,
                      datecol = 'datetime',
                      agg = 'annual',
                      minNumObs = 100,
                      minNumUncen = 50
                     )

                    wrtds_out <- flux_annual_wrtds %>%
                         filter(wy %in% good_years) %>%
                         rename(val = flux) %>%
                         mutate(site_code = site_code,
                              var = solutes[j],
                              method = 'wrtds',
                              ms_recommended = 0)

                        return(wrtds_out)
                  },
                  error = function(e) {
                    writeLines(paste('\nWRTDS run failed for \n     site', site_code,
                                     '\n     variable', target_solute, '\n WRTDS TRYING AGAIN'))
                    tryCatch(
                      expr = {
                        flux_annual_wrtds <- calculate_wrtds(
                               chem_df = con_full,
                               q_df = q_df,
                               ws_size = area,
                               lat = lat,
                               long = long,
                               datecol = 'datetime',
                               agg = 'annual',
                               minNumObs = 100,
                               minNumUncen = 50
                        )

                       wrtds_out <- flux_annual_wrtds %>%
                         filter(wy %in% good_years) %>%
                         rename(val = flux) %>%
                         mutate(site_code = site_code,
                              var = solutes[j],
                              method = 'wrtds',
                              ms_recommended = 0)

                        return(wrtds_out)
                      },
                      error = function(e) {
                        print("WRTDS failed, setting to NA")
                        flux_annual_wrtds <- NA
                      }
                    )
                  }
                ) # end wrtds

              }

              # if not wrtdsk
              for(k in 1:length(good_years)){

               writeLines(paste("site:", site_code,
                              'year:', good_years[k]))

               target_year <- as.numeric(as.character(good_years[k]))

               # calculate flag ratios to carry forward
               flag_df <- carry_flags(raw_q_df = raw_data_q,
                                      raw_con_df = raw_data_con_in,
                                      target_year = target_year,
                                      target_solute = target_solute,
                                      period = 'annual')

               raw_data_target_year <- raw_data_full %>%
                   mutate(wy = as.numeric(as.character(wy))) %>%
                   filter(wy == target_year)

               q_target_year <- raw_data_target_year %>%
                   select(site_code, datetime, q_lps, wy)%>%
                   na.omit()

               con_target_year <- raw_data_target_year %>%
                   select(site_code, datetime, con, wy) %>%
                   na.omit()

               ### calculate annual flux ######
               chem_df_errors <- con_target_year
               q_df_errors <- q_target_year

               ### save and then remove errors attribute for calcs
               chem_df <- errors::drop_errors(chem_df_errors)
               q_df <- errors::drop_errors(q_df_errors)

              if(method == 'average') {
               #### calculate average ####
                flux_annual <- raw_data_target_year %>%
                   group_by(wy) %>%
                   summarize(q_lps = mean(q_lps, na.rm = TRUE),
                             con = mean(con, na.rm = TRUE)) %>%
                   # multiply by seconds in a year, and divide my mg to kg conversion (1M)
                   mutate(flux = con*q_lps*3.154e+7*(1/area)*1e-6) %>%
                 pull(flux)
              } else if(method == 'pw') {
                 #### calculate period weighted #####
                 flux_annual <- calculate_pw(chem_df, q_df, datecol = 'datetime')
              } else if (method =='beale') {
                 #### calculate beale ######
                 flux_annual <- calculate_beale(chem_df, q_df, datecol = 'datetime')

              } else if (method == 'rating') {
                 #### calculate rating #####
                 flux_annual <- calculate_rating(chem_df, q_df, datecol = 'datetime')
              } else if (method == 'composite') {
                  #### calculate composite ######
                  rating_filled_df <- generate_residual_corrected_con(chem_df = chem_df,
                                                                      q_df = q_df,
                                                                      datecol = 'datetime',
                                                                      sitecol = 'site_code')
                  # calculate annual flux from composite
                  flux_annual_comp <- calculate_composite_from_rating_filled_df(rating_filled_df)
                  flux_annual <- flux_annual_comp$flux[1]
              }

                 #### congeal fluxes ####
                 target_year_out <- tibble(wy = as.character(target_year),
                                           val = flux_annual,
                                           site_code = !!site_code,
                                           var = !!target_solute,
                              method = !!method)
                 out_frame <- bind_rows(out_frame, target_year_out)
              } # end year loop

            } else {

        # back to simple flux
        chem_is_highres <- Mode(diff(as.numeric(chem_chunk$datetime))) <= 15 * 60
        if(is.na(chem_is_highres)) { chem_is_highres <- FALSE}

                            #if both chem and flow data are low resolution (grab samples),
                            #   let approxjoin_datetime match up samples with a 12-hour gap. otherwise the
                            #   gap should be 7.5 mins so that there isn't enormous duplication of
                            #   timestamps where multiple high-res values can be snapped to the
                            #   same low-res value
                            if(! chem_is_highres && ! flow_is_highres) {
                                join_distance <- c('12:00:00')#, '%H:%M:%S')
                            } else {
                                join_distance <- c('7:30')#, '%M:%S')
                            }

                            chem_split[[i]] <- approxjoin_datetime(x = chem_chunk,
                                                                   y = site_q,
                                                                   rollmax = join_distance,
                                                                   keep_datetimes_from = 'x')

                          if(q_type == 'discharge'){
                            chem_split[[i]] <- chem_split[[i]] %>%
                              mutate(site_code = site_code_x,
                                     var = var_x,
                                     # kg/interval = mg/L *  L/s  * q_interval / 1e6
                                     val = val_x * val_y * errors::as.errors(q_interval) / errors::as.errors(1e6),
                                        ms_status = numeric_any_v(ms_status_x, ms_status_y),
                                        ms_interp = numeric_any_v(ms_interp_x, ms_interp_y)) %>%
                                 select(-starts_with(c('site_code_', 'var_', 'val_',
                                                       'ms_status_', 'ms_interp_'))) %>%
                                 filter(! is.na(val)) %>% #should be redundant
                               arrange(datetime) %>%
                               ms_scale_flux_by_area(site_info)
                            } else {
                                chem_split[[i]] <- chem_split[[i]] %>%
                                    mutate(site_code = site_code_x,
                                           var = var_x,
                                           # kg/interval/ha = mg/L *  mm/interval * ha/100
                                           val = val_x * val_y / errors::as.errors(100),
                                           ms_status = numeric_any_v(ms_status_x, ms_status_y),
                                           ms_interp = numeric_any_v(ms_interp_x, ms_interp_y)) %>%
                                    select(-starts_with(c('site_code_', 'var_', 'val_',
                                                          'ms_status_', 'ms_interp_'))) %>%
                                    filter(! is.na(val)) %>% #should be redundant
                                    arrange(datetime)
                            }

                        flux <- chem_split %>%
                            purrr::reduce(bind_rows) %>%
                            arrange(site_code, var, datetime)

                        all_sites_flux <- rbind(all_sites_flux, flux)


                    if(nrow(all_sites_flux) == 0) { return(NULL) }

                    all_sites_flux$val_err <- errors::errors(all_sites_flux$val)
                    all_sites_flux$val <- errors::drop_errors(all_sites_flux$val)

                    return(all_sites_flux)
                }
            }
    }

    return(out_frame)
}

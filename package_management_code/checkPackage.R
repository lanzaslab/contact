devtools::document() #add/update man files

devtools::check(remote = TRUE, run_dont_test = TRUE, vignettes = FALSE) #check on local system

devtools::check_win_release() #check on windows server running the current version of R

devtools::check_win_devel() #check on windows server running the development version of R

devtools::revdep() #check for any downstream dependencies

devtools::build() #build tarball


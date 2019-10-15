
#tell the build to ignore specific folders
usethis::use_build_ignore("data-raw")
usethis::use_build_ignore("package_management_code") #folder containing package-building directions.
usethis::use_build_ignore("inProgress_functions") #folder containing functions that we are currently building and will be added to later versions of the package.
#usethis::use_build_ignore("vignettes")

devtools::build_vignettes(pkg = ".", dependencies = "VignetteBuilder",
                clean = TRUE, upgrade = "never", quiet = FALSE, install = FALSE,
                keep_md = TRUE) #Create vignette

usethis::use_readme_rmd() #create README.Rmd file

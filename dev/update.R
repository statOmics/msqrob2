## ********************************************
## Update your package code before a git commit
## ********************************************
## This template was modified from https://lcolladotor.github.io/biocthis/

## Automatically re-style the code in your package to a Bioconductor-friendly
## format
styler::style_pkg(transformers = biocthis::bioc_style())
styler::style_dir(usethis::proj_path("dev"), transformers = biocthis::bioc_style())
styler::style_dir(
    usethis::proj_path("vignettes"),
    transformers = biocthis::bioc_style(),
    filetype = "Rmd"
)

## Re-make the documentation files
devtools::document()

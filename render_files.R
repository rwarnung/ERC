## more style: https://jozef.io/r913-spin-with-style/
## themese: https://www.datadreaming.org/post/r-markdown-theme-gallery/
knitr::spin("./ERC_optim.R", knit = FALSE)
#rmarkdown::render("./R_example.R")

rmarkdown::render(
  "./ERC_optim.R", 
  output_format = rmarkdown::html_document(
    theme = "united"
    #,
    #mathjax = NULL,
    #highlight = NULL
  )
)
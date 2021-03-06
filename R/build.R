# passed from blogdown::build_site()
local <- commandArgs(TRUE)[1] == "TRUE"

knitr::opts_knit$set(
  base.dir = normalizePath("static/", mustWork = TRUE),
  base.url = "/"
)

rmds <- list.files("content", "[\\.]Rmd$", recursive = TRUE, full.names = TRUE)

for (rmd in rmds) {
  wo_ext <- tools::file_path_sans_ext(rmd)
  html <- glue::glue("{wo_ext}.html")
  
  if (file.exists(html) && utils::file_test("-ot", rmd, html)) {
    message(glue::glue("skip {rmd}"))
    next
  }
  
  knitr::opts_chunk$set(
    fig.path = glue::glue("post/{basename(wo_ext)}_files/figure-html/")
  )
  
  set.seed(1984)
  # knitr::knit(input = rmd, output = html, encoding = "UTF-8")
  blogdown:::build_rmds(rmd)
}

blogdown::hugo_build(local = local)
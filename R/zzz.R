.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Thank you for using microbAIDeR! For detailed documentation, visit https://github.com/FabbriniMarco/microbAIDeR/wiki.\n\n",
                        paste(capture.output(citation("microbAIDeR")), collapse = "\n"))
}

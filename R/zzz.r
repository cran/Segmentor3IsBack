.onAttach <- function(lib, pkg) {
    library.dynam("Segmentor3IsBack", pkg, lib)
    packageStartupMessage("Segmentor3IsBack Loaded \n")
   
    
}



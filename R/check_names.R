# Check plant species names for typos and synonyms
#' \code{check_names} to check plant species names in Wisconsin for typos and synonyms using UWSP herbarium as baseline.
#' 
#' @author Daijiang Li
#' 
#' @param sp.to.check A vector of plant species to check.
#' @param sp.dataset species names and their synonyms from UWSP herbarium website.
#' @examples
#' check_names(sp.to.check = c("Acer rubrum", "Quercus abla"))
#' @return
#' \itemize{
#'   \item sp.origin the original species names.
#'   \item sp.spell the closest correct species names after check for typos, if any.
#'   \item sp.final the updated species names after check for synonyms, if any.
#'   \item spell.correct logical, are the species names spelled correctly?
#'   \item syc.correct logical, are the species names updated?
#' }
#' @export
check_names <- function(sp.to.check, sp.dataset = sp_syn_uwsp){
  sp.dataset.list = unique(na.omit(sp.dataset$syn))
  #most close spelling in the dataset
  sp.spell.checked = sapply(sp.to.check, function(x){
    # igonre xxx.spp xxx sp
    if(any(str_detect(pattern = ".*sp[0-9]?$|.*spp[0-9]?$|.*sp\\.$", string = x))){
      x
    } else{
      sp.dataset.list[which.min(adist(x, sp.dataset.list))]}
  })
  #sp names used as main in webpages of UWSP herbarium
  sp.syc = sapply(sp.spell.checked, function(x){
    if(any(str_detect(pattern = ".*sp[0-9]?$|.*spp[0-9]?$|.*sp\\.$", string = x))){
      x
    } else{
      y = unique(na.omit(sp.dataset$sp[sp.dataset$syn == x]))
      if(length(y > 1)){
        y = y[which.min(adist(x, y))]
      }
      y}
  })
  data.frame(sp.origin = sp.to.check, sp.spell = sp.spell.checked,
             sp.final = sp.syc, 
             spell.correct = sp.to.check == sp.spell.checked,
             syc.correct = sp.spell.checked == sp.syc)
}

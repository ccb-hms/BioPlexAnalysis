############################################################
# 
# author: Ludwig Geistlinger
# date: 2023-01-13 11:19:25
# 
# descr: functionality for working with protein structure
#        data from PDB
# 
############################################################

#' View a protein structure from PDB
#' 
#' @description This function spins up a 3D viewer of a protein structure
#' @param pdb.id character. A PDB ID in 4-letter PDB code or 6-letter
#' PDB-ID_Chain-ID code.
#' @return An html widget that can be explored in the browser. 
#' @examples
#'  viewPDB("6nmi")
#'
#' @importFrom r3dmol "%>%"
#' @export
viewPDB <- function(pdb.id)
{
    # obtain pdb file
    stopifnot(is.character(pdb.id) && length(pdb.id) == 1)
    pdb.file <- bio3d::get.pdb(pdb.id)
    pdb <- bio3d::read.pdb(pdb.file)
    file.remove(pdb.file)

    # color the chains
    chains <- unlist(pdb$remark$biomat$chain)
    nr.chains <- length(chains)
    chain.colors <- RColorBrewer::brewer.pal(nr.chains, "Set1")
    chain.colors <- tolower(chain.colors)

    viewer <- r3dmol::r3dmol(id = "", elementId = "demo") %>%
        # Add model to scene
        r3dmol::m_add_model(data = r3dmol::m_bio3d(pdb), format = "pdb") %>%
        # Zoom to encompass the whole scene
        r3dmol::m_zoom_to() %>%
        # Set style of structures
        r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "#00cc96"))
        
    for(i in seq_len(nr.chains))
        viewer <- r3dmol::m_set_style(viewer,
                              sel = r3dmol::m_sel(chain = LETTERS[i]),
                              style = r3dmol::m_style_cartoon(color = chain.colors[i])) 

    viewer %>% r3dmol::m_rotate(angle = 90, axis = "y") %>% r3dmol::m_spin()
}

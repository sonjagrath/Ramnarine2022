#' Associates modules with sample information
#'
#' @description We want to know which modules are associated with interesting sample characteristics.
#' This is achieved by correlating the sample property, i.e. stress condition or control,
#' to the module eigengene
#'
#' @param inf
#' @param MEs
#' @param type
#' @param cortype
#' @param lab
#'
#' @return
#' @export
#'
#' @examples
Module2SampleFeat <- function(inf, MEs, type = "categorical" , cortype = 'pearson', lab = "" ){
    require(broom)
    df<- inf %>% dplyr::bind_cols(MEs) %>%
        tidyr::gather("VarN","varID", 1:dim(inf)[2]) %>%
        tidyr::gather("module","eigenvalue", -VarN, -varID) %>%
        dplyr::group_by(module,VarN)

    if(  type == "categorical"){
        # for categorical variables do lm (anova)
        df<- df %>% dplyr::do( tidy( lm(scale(eigenvalue) ~ varID,.))) %>%
            dplyr::filter(term != "(Intercept)") %>%
            dplyr::mutate( mcol = sub( "ME", "", module),
                           infVar = paste( VarN, sub("varID","",term), sep = "\n"  ))
    }else{
        # for continuous calculate Pearsons correlation
        df<- df %>% dplyr::do( tidy( cor.test(.$eigenvalue,.$varID,method = cortype))) %>%
            dplyr::mutate(mcol = sub( "ME", "", module),
                          infVar = VarN)
    }

    df$mcol<-factor(df$mcol)
    mc<-ggplot(df,aes(x="A",y=module,fill=module))+
        geom_tile()+
        scale_fill_manual(values = levels( df$mcol))+
        theme(legend.position ="None",
              axis.title = element_blank(),
              axis.text.x=element_blank())

    tp<- ggplot(df,aes( x=infVar, y=module, fill=estimate,label=format(p.value,scientific=T,digits=2) )) +
        geom_tile()+
        geom_text()+
        scale_fill_gradient2()+
        theme(#axis.title = element_blank(),
            legend.position = "bottom",
            legend.key.width=unit(1.7,"cm"),
            axis.text.y=element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_line(size=NA),
            axis.text.x = element_text(angle =45 ,hjust = 1))

    cowplot::plot_grid(mc,tp, nrow = 1, rel_widths = c(0.3,0.7),align = "h",axis="lb",labels =lab )
}

#' Function to get all genes from specific module
#'
#' @param color
#'
#' @return
#' @export
#'
#' @examples
get_kim_color <- function(color){
    kim_color <- kim[moduleColors == color, ]
    kim_color <- kim_color[order(kim_color$kWithin, decreasing = TRUE), ]
    kim_color <- data.frame(ENSEMBL = rownames(kim_color), kim_color[ , 1:2],
                            stringsAsFactors = FALSE)
}


#' Get specific number of hubs
#'
#' @param color
#' @param n
#'
#' @return
#' @export
#'
#' @examples
get_hubs <- function(color, n){
    km <- AnnotationDbi::select(org.Dm.eg.db,
                                keys = data$ENSEMBL,
                                keytype = "ENSEMBL",
                                columns = c("ENSEMBL", "SYMBOL")) %>%
        right_join(get_kim_color(color))
    km[1:n,]
}


#' Combine DESeq2 results with hubs
#'
#' @param results
#' @param color
#' @param n
#'
#' @return
#' @export
#'
#' @examples
are_hubs_DE <- function(results, color, n){
    results[get_kim_color(color)$ENSEMBL[1:n],]
}


#' Annotate a specific module from WGCNA
#'
#' @param color
#'
#' @return
#' @export
#'
#' @examples
annotate_module <- function(color){
    g4reac %>%
        dplyr::filter(module == color) %>%
        na.omit() %>%
        dplyr::select(ENTREZID) %>%
        unique()
}


#' Enrich specific module from WGCNA
#'
#' @param color
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
module_enrich_pathway <- function(color, organism = "fly"){
    enrichPathway(gene = annotate_module(color)$ENTREZID,
                  organism = organism,
                  readable = TRUE,
                  pvalueCutoff = 0.05)
}


#' Complete visualization for specific module
#'
#' @param color
#' @param nr_categories
#' @param type
#' @param path
#'
#' @return
#' @export
#'
#' @examples
complete_visualization <- function(color,
                                   nr_categories,
                                   type = "dotplot",
                                   path = getwd()){
    x <- module_enrich_pathway(color, organism = "fly")

    if (type == "dotplot") {
        p1 <- dotplot(x, showCategory = nr_categories) +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 40))
        ggsave(filename = paste(path, "/", color, "/", color,"_dotplot",".png",
                                sep=""),
               width = 10,
               height = 10)
    }

    if (type == "emapplot") {
        p2 <- emapplot(x, layout="kk", showCategory = nr_categories)
        ggsave(filename = paste(path, "/", color, "/", color,"_emapplot",".png",
                                sep=""),
               width = 10,
               height = 10)
    }

    if (type == "cnetplot") {
        p3 <- cnetplot(x,
                       foldChange = geneList,
                       layout = "kk",
                       # foldChange = geneList[geneList>=2],
                       showCategory = nr_categories,
                       node_label = "all")
        ggsave(filename = paste(path, "/", color, "/", color,"_cnetplot",".png",
                                sep=""),
               width = 10,
               height = 10)
    }
    # TODO:
    # Add possibility to give several plot variants in 'type' argument
    ## all_plots <- list(p1, p2, p3)
    ## return(all_plots)
}


#' Creates subdirectory for given color
#'
#' @param path
#' @param color
#'
#' @return
#' @export
#'
#' @examples
create_color_folder <- function(path, color){
    dir.create(path = paste(path, "/", color, sep = ""))
}


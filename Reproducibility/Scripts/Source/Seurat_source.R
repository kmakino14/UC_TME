suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(viridisLite))
suppressMessages(library(BuenColors))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(harmony))
suppressMessages(library(pheatmap))
#suppressMessages(library(Nebulosa))
suppressMessages(library(mclust))
suppressMessages(library(future))
suppressMessages(library(ggrastr))
suppressMessages(library(ggrepel))

options(stringsAsFactors=F)

######################################################
################# general function          ##########
######################################################

# color_Palette = readRDS("~/data_home/colorPalette.rds")

take_factor =  function(list,order=order,sep=sep){sapply(strsplit(list,sep),function(x){if(length(order)==1){x[order]}else{paste(x[order],collapse=sep)}})}

write.table_FT_2 = function(dataframe,filename="tmp.txt"){dirname(filename)%>%dir.create_p();write.table(dataframe,filename,row.names=F, col.names=T, sep="\t", append=F, quote=F)}
write.table_n_2  = function(dataframe,name_of_row="rowname",filename="tmp.txt"){dirname(filename)%>%dir.create_p();dataframe = dataframe %>% rownames_to_column(name_of_row);write.table(dataframe,filename,row.names=F, col.names=T, sep="\t", append=F, quote=F)}

read.table_FT = function(filename){read.table(filename,header=T,row.names=NULL,stringsAsFactors=F)}

fread_n   = function(filename){data<-fread(filename);data<-as.data.frame(data);rownames(data)<-data[,1];data<-data[,-1,drop=F];data}
fread_FT  = function(filename){data<-fread(filename);data<-as.data.frame(data);data}

pdf_2  = function(filename,h=5,w=5){dirname(filename)%>%dir.create_p();pdf(filename,h=h,w=w)}
png_2  = function(filename,h=480,w=480,pointsize=20,res=300){dirname(filename)%>%dir.create_p();png(filename,height=h,width=w,pointsize=pointsize,res=res)}
plot_2 = function(filename,nrow=1,ncol=1,pointsize=20){if(grepl("pdf$",filename)){pdf_2(filename,h=5*nrow,w=5*ncol)}
                                                       if(grepl("png$",filename)){png_2(filename,h=480*nrow,w=480*ncol,pointsize=pointsize)}
                                                      }
## h=960, w = 960

"%ni%" <- Negate("%in%")

dir.create_p = function(dirname){dir.create(dirname, showWarnings = FALSE, recursive = TRUE)}

hclust_order = function(DATA,METHOD="ward.D2"){
  suppressMessages(library(ggdendro))
  rd     = dist(DATA)
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  col_order = as.character(ddata$labels$label)
  col_order
}

##########################
plot_theme    = function(size=10,theme="classic",legend=TRUE,x_angle=90){
                                     if(theme == "classic"){p = theme_classic()}
                                     if(theme == "void"){   p = theme_void()}
                                     if(x_angle==90){p =  p + theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = x_angle , hjust = 1 ,  vjust = 0.5 ))}else{
                                      if(x_angle==0){p =  p + theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = x_angle , hjust = 0.5 ,  vjust = 0.5 ))}else{
                                                     p =  p + theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = x_angle , hjust = 1 , vjust = 1 ))
                                                    }}
                                     p = p + theme(axis.text.y  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
                                                   axis.title.x = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
                                                   axis.title.y = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 90),
                                                   plot.title   = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", hjust = 0.5))
                                     if(legend){ p = p + theme(legend.title = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
                                                               legend.text  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))}else{
                                                 p = p + theme(legend.position  = "none")
                                               }
                                     return(p)
                                    }

#### plot_classic    = function( size=10){ p = theme_classic()+
####                                            theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 60 , hjust = 1 , vjust = 1 ),
####                                                  axis.text.y  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                  axis.title.x = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                  axis.title.y = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 90),
####                                                  plot.title   = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                  legend.title = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                  legend.text  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))
####                                      return(p)
####                                     }
#### 
#### plot_classic_v = function( size=10){ p = theme_classic()+
####                                           theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                 axis.text.y  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                 axis.title.x = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                 axis.title.y = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 90),
####                                                 plot.title   = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                 legend.title = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                 legend.text  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))
####                                      return(p)
####                                     }
#### 
#### plot_classic_wo_leg = function(size=10){ p = theme_classic()+
####                                               theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                     axis.text.y  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                     axis.title.x = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                     axis.title.y = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 90),
####                                                     plot.title   = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                                     legend.position  = "none" )
####                                          return(p)
####                                        }
#### 
#### plot_void = function( plot = p, size=13){ plot + theme_void()+
####                                       theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 60 , hjust = 1 , vjust = 1 ),
####                                             axis.text.y  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                             axis.title.x = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                             axis.title.y = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 90),
####                                             plot.title   = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                             legend.title = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
####                                             legend.text  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))
####                                  }

data_max_min = function(data = data , max = NA , min = NA){data[data > max] = max
                                                           data[data < min] = min
                                                           data}

heatmap_col = c("#6594C4", "#77A6CE", "#8BB9D7", "#9EC7DF", "#B2D4E7", "#C6E2EE", "#DAEEF5",
                "#E5F5EC", "#EDF8DE", "#F4FBD0", "#FDFDC1", "#FEF8B5", "#FEF0A9", "#FEE89D",
                "#FEE091", "#FDCE83", "#FDB976", "#FCA468", "#FC8F5A", "#F4794E")

heatmap_gap = function(target = target, level = level){gaps_res=c()
                                                       for(tmp_name in level){
                                                         tmp_num = grep(tmp_name,target) %>% max(.)
                                                         gaps_res = c(gaps_res,tmp_num)
                                                       }
                                                       gaps_res
                                                      }

correlation_col = c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC", "#053061") %>%
                   rev(.) %>% colorRampPalette(.)

bed_arrange = function(data=data){data$Chr   = factor(data$Chr,levels=paste0("chr",c(1:22,"X","Y")))
                                  data$Start = as.numeric(data$Start)
                                  res = data %>% dplyr::arrange(Chr,Start)
                                  res}

rowname_clear = function(data){rownames(data)=NULL;data}

######################################################
################# ASAP QC plot         ###############
######################################################

ASAP_QC_plot = function(ASAP   = ASAP,
                        OUTPUT = paste0("plot_QC_plot/","tmp","_","ASAP","_","HTOcondition","_QCplot.png"),
                        group.by = "HTO_final_2",
                        title = "",
                        pt.size= 0,
                        ncol = 3,
                        w_enlarge=1,
                        pct_reads_in_peaks_tmp   = 30,
                        peak_region_fragment_tmp = 50000,
                        passed_filters_tmp       = 1000,
                        passed_filters_tmp2      = 100000,
                        TSS.enrichment_tmp       = 2,
                        nucleosome_signal_tmp    = 2.5,
                        nCount_ADT_tmp           = 500,
                        nCount_ADT_tmp2         = 30000
                        ){
                         # Visualize QC metrics as a violin plot
                         p=list()
                         p[[1]]  = VlnPlot(ASAP, c("pct_reads_in_peaks"),   pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + 
                                    geom_hline(yintercept=pct_reads_in_peaks_tmp ,  col="black")
                         p[[2]]  = VlnPlot(ASAP, c("peak_region_fragments"),pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() + 
                                    geom_hline(yintercept=peak_region_fragment_tmp,col="black")
                         p[[3]]  = VlnPlot(ASAP, c("passed_filters"),       pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=passed_filters_tmp, col="black") +
                                    geom_hline(yintercept=passed_filters_tmp2,col="black")
                         p[[4]]  = VlnPlot(ASAP, c("TSS.enrichment"),       pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() +
                                    geom_hline(yintercept=TSS.enrichment_tmp,    col="black")
                         p[[5]]  = VlnPlot(ASAP, c("nucleosome_signal"),    pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() +
                                    geom_hline(yintercept=nucleosome_signal_tmp,  col="black")
                         p[[6]]  = VlnPlot(ASAP, c("nCount_ADT"),           pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=nCount_ADT_tmp,  col="black") +
                                    geom_hline(yintercept=nCount_ADT_tmp2,col="black")
                         
                         p_sum = patchwork::wrap_plots(p, ncol = ncol)+
                                  plot_annotation(paste0(title," QC ,",ncol(ASAP),"cells"))

                         nrow = 6/ncol %>% ceiling(.)
                         
                         if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow,w=5*ncol*w_enlarge)}
                         if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow,w=480*ncol*w_enlarge)}
                          plot(p_sum)
                         dev.off()
                        }

######################################################
################# CITE QC plot         ###############
######################################################

CITE_QC_plot = function(Seurat = CITE,
                        OUTPUT = paste0("plot_QC_plot/","tmp","_","CITE","_","HTOcondition","_QCplot.png"),
                        group.by = "HTO_final_2",
                        title = "",
                        pt.size= 0,
                        ncol = 4,
                        w_enlarge=1,
                        percent.mt_tmp   = 20,
                        nFeature_RNA_tmp  = 1000,
                        nFeature_RNA_tmp2 = 7000,
                        nCount_RNA_tmp    = 1000,
                        nCount_RNA_tmp2   = 30000,
                        nCount_ADT_tmp    = 200,
                        nCount_ADT_tmp2  = 10000
                        ){
                         # Visualize QC metrics as a violin plot
                         p=list()
                         p[[1]]  = VlnPlot(Seurat, c("percent.mt"),   pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + 
                                    geom_hline(yintercept=percent.mt_tmp ,  col="black")
                         p[[2]]  = VlnPlot(Seurat, c("nFeature_RNA"), pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() + 
                                    geom_hline(yintercept=nFeature_RNA_tmp,col="black")+
                                    geom_hline(yintercept=nFeature_RNA_tmp2,col="black")
                         p[[3]]  = VlnPlot(Seurat, c("nCount_RNA"),   pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=nCount_RNA_tmp, col="black") +
                                    geom_hline(yintercept=nCount_RNA_tmp2,col="black")
                         p[[4]]  = VlnPlot(Seurat, c("nCount_ADT"),   pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=nCount_ADT_tmp,  col="black") +
                                    geom_hline(yintercept=nCount_ADT_tmp2,col="black")
                         
                         p_sum = patchwork::wrap_plots(p, ncol = ncol)+
                                  plot_annotation(paste0(title," QC ,",ncol(Seurat),"cells"))

                         nrow = 4/ncol %>% ceiling(.)
                         
                         if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow,w=5*ncol*w_enlarge)}
                         if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow,w=480*ncol*w_enlarge)}
                          plot(p_sum)
                         dev.off()
                        }

######################################################
################# icCITE QC plot       ###############
######################################################

icCITE_QC_plot = function(Seurat = icCITE,
                          OUTPUT = paste0("plot_QC_plot/","tmp","_","icCITE","_","HTOcondition","_QCplot.png"),
                          group.by = "HTO_final_2",
                          title = "",
                          pt.size= 0,
                          ncol = 4,
                          w_enlarge=1,
                          percent.mt_tmp   = 20,
                          nFeature_RNA_tmp  = 1000,
                          nFeature_RNA_tmp2 = 7000,
                          nCount_RNA_tmp    = 1000,
                          nCount_RNA_tmp2   = 30000,
                          nCount_ADT_tmp    = 200,
                          nCount_ADT_tmp2  = 10000,
                          nCount_TSB_tmp    = 200,
                          nCount_TSB_tmp2  = 10000
                        ){
                         # Visualize QC metrics as a violin plot
                         p=list()
                         p[[1]]  = VlnPlot(Seurat, c("percent.mt"),   pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + 
                                    geom_hline(yintercept=percent.mt_tmp ,  col="black")
                         p[[2]]  = VlnPlot(Seurat, c("nFeature_RNA"), pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() + 
                                    geom_hline(yintercept=nFeature_RNA_tmp,col="black")+
                                    geom_hline(yintercept=nFeature_RNA_tmp2,col="black")
                         p[[3]]  = VlnPlot(Seurat, c("nCount_RNA"),   pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=nCount_RNA_tmp, col="black") +
                                    geom_hline(yintercept=nCount_RNA_tmp2,col="black")
                         p[[4]]  = VlnPlot(Seurat, c("nCount_ADT"),   pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=nCount_ADT_tmp,  col="black") +
                                    geom_hline(yintercept=nCount_ADT_tmp2,col="black")
                         p[[5]]  = VlnPlot(Seurat, c("nCount_TSB"),   pt.size = pt.size, group.by = group.by) + 
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept=nCount_TSB_tmp,  col="black") +
                                    geom_hline(yintercept=nCount_TSB_tmp2,col="black")

                         p_sum = patchwork::wrap_plots(p, ncol = ncol)+
                                  plot_annotation(paste0(title," QC ,",ncol(Seurat),"cells"))

                         nrow = 5/ncol %>% ceiling(.)
                         
                         if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow,w=5*ncol*w_enlarge)}
                         if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow,w=480*ncol*w_enlarge)}
                          plot(p_sum)
                         dev.off()
                        }

######################################################
################# DOGMA QC plot         ###############
######################################################

DOGMA_QC_plot = function(Seurat   = DOGMA,
                         OUTPUT = paste0("plot_QC_plot/","tmp","_","CITE","_","HTOcondition","_QCplot.png"),
                         group.by = "HTO_final_2",
                         title = "",
                         pt.size= 0,
                         ncol = 9,
                         w_enlarge=1,
                         pct_reads_in_peaks_tmp   = 20,
                         peak_region_fragment_tmp = 35000,
                         passed_filters_tmp       = 500,
                         passed_filters_tmp2      = 60000,
                         TSS.enrichment_tmp       = 2,
                         nucleosome_signal_tmp    = 2.5,
                         percent.mt_tmp           = 30,
                         nFeature_RNA_tmp         = 500,
                         nFeature_RNA_tmp2        = 6000,
                         nCount_RNA_tmp           = 500,
                         nCount_RNA_tmp2          = 20000,
                         nCount_ADT_tmp           = 1000,
                         nCount_ADT_tmp2          = 10000
                        ){
                         # Visualize QC metrics as a violin plot
                         p=list()
                         p[[1]]  = VlnPlot(Seurat, c("pct_reads_in_peaks"),         pt.size = pt.size, group.by = group.by) +
                                    NoLegend()                   +
                                    geom_hline(yintercept = pct_reads_in_peaks_tmp ,  col="black")
                         p[[2]]  = VlnPlot(Seurat, c("atac_peak_region_fragments"), pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept = peak_region_fragment_tmp ,col="black")
                         p[[3]]  = VlnPlot(Seurat, c("atac_fragments"),             pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept = passed_filters_tmp, col="black") + 
                                    geom_hline(yintercept = passed_filters_tmp2,col="black")
                         p[[4]]  = VlnPlot(Seurat, c("TSS.enrichment"),             pt.size = pt.size, group.by = group.by) +
                                    NoLegend()                   +
                                    geom_hline(yintercept = TSS.enrichment_tmp,    col="black")
                         p[[5]]  = VlnPlot(Seurat, c("nucleosome_signal"),          pt.size = pt.size, group.by = group.by) +
                                    NoLegend()                   +
                                    geom_hline(yintercept = nucleosome_signal_tmp,  col="black")
                         p[[6]]  = VlnPlot(Seurat, c("percent.mt"),                 pt.size = pt.size, group.by = group.by) +
                                    NoLegend()                   +
                                    geom_hline(yintercept = percent.mt_tmp ,   col="black")
                         p[[7]]  = VlnPlot(Seurat, c("nFeature_RNA"),               pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept = nFeature_RNA_tmp,   col="black") +
                                    geom_hline(yintercept = nFeature_RNA_tmp2, col="black")
                         p[[8]]  = VlnPlot(Seurat, c("nCount_RNA"),                 pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept = nCount_RNA_tmp,  col="black") +
                                    geom_hline(yintercept = nCount_RNA_tmp2, col="black")
                         p[[9]]  = VlnPlot(Seurat, c("nCount_ADT"),                 pt.size = pt.size, group.by = group.by) +
                                    NoLegend() + scale_y_log10() +
                                    geom_hline(yintercept = nCount_ADT_tmp,  col="black") +
                                    geom_hline(yintercept = nCount_ADT_tmp2, col="black")
                         
                         p_sum = patchwork::wrap_plots(p, ncol = ncol)+
                                  plot_annotation(paste0(title," QC ,",ncol(Seurat),"cells"))

                         nrow = 9/ncol %>% ceiling(.)
                         
                         if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow,w=5*ncol*w_enlarge)}
                         if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow,w=480*ncol*w_enlarge)}
                          plot(p_sum)
                         dev.off()
                         }


######################################################
################# plot of density plot ###############
######################################################
plot_density = function(Seurat = data,
                        UMAP   = "wnn.RNA.ADT.umap",
                        ident  = "HTO_final_2",
                        target_list = HTO_list_2$HTO,
                        OUTPUT = OUTPUT,
                        name   = name,
                        ncol = NA,
                        plot_thresh=1.1,
                        w_enlarge = 1,
                        h_enlarge = 1,
                        save=TRUE
                        ){
                   suppressMessages(library(RColorBrewer))
                   suppressMessages(library(tidyverse))
                   suppressMessages(library(ggrastr))

                   Idents(Seurat)      = ident
                   
                   UMAPdata = Embeddings(Seurat, reduction = UMAP) %>% as.data.frame() %>% rownames_to_column("CellBarcode") %>%
                               left_join(.,Seurat@meta.data %>% as.data.frame() %>% rownames_to_column("CellBarcode"),by="CellBarcode")
                   colnames(UMAPdata)[2:3] = c("UMAP_1","UMAP_2")
                   
                   col_list = colorRampPalette(jdb_palette("brewer_heat")[3:9])
                   col=col_list(50)
                   
                   x_min = min(UMAPdata$UMAP_1) * plot_thresh
                   x_max = max(UMAPdata$UMAP_1) * plot_thresh
                   y_min = min(UMAPdata$UMAP_2) * plot_thresh
                   y_max = max(UMAPdata$UMAP_2) * plot_thresh
                   
                   ##########################
                   
                   p=list()
                   for(iii in 1:length(target_list)){
                       tmp_target   = target_list[iii]
                       p[[tmp_target]] = ggplot(UMAPdata, aes(x=UMAP_1,y=UMAP_2))+
                                   stat_density_2d(geom = "polygon",alpha=0.3,fill="#FEE8C8",color="gray60",breaks=0.002)+
                                   geom_point_rast(size=.1,color="gray")+
                                   stat_density_2d(data =UMAPdata[UMAPdata[,ident]==tmp_target,],
                                                   aes(x=UMAP_1,y=UMAP_2,fill = after_stat(level)), geom = "polygon",bins=10)+
                                   scale_fill_gradientn(colours=col)+
                                   xlim(c(x_min,x_max))+
                                   ylim(c(y_min,y_max))+
                                   plot_theme(theme="void",legend=FALSE,size=0)+
                                   theme(plot.title   = element_text(colour = "black", size = 10, face = "bold", family = "Helvetica", hjust = 0.5))+
                                   labs(title=tmp_target)
                       print(tmp_target)
                   }
                   
                   if(save){
                      if(is.na(ncol)){ncol_tmp = sqrt(length(target_list) * 2 ) %>% ceiling(.)}else{ncol_tmp=ncol}
                      nrow_tmp = (length(target_list)/ncol_tmp) %>% ceiling(.)
                      
                      p_sum = patchwork::wrap_plots(p, ncol = ncol_tmp)+
                               plot_annotation(paste0(SAMPLE_TMP," : ",name," : ","HTO"))
                      
                      ## p[[1]] = AugmentPlot(plot=p[[1]],width=4,height=4,dpi=500)
                      
                      if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow_tmp*h_enlarge,w=5*ncol_tmp*w_enlarge)}
                      if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow_tmp*h_enlarge,w=480*ncol_tmp*w_enlarge)}
                       plot(p_sum)
                      dev.off()
                   }else{p}
                   }

######################################################
################# plot of density Nebulosa ###############
######################################################
plot_density_2 = function(Seurat = data,
                          UMAP   = "wnn.RNA.ADT.umap",
                          target = "HTO_final",
                          target_list = NA , 
                          OUTPUT = OUTPUT,
                          SAMPLE_TMP = SAMPLE_TMP,
                          name   = "name",
                          palette = "viridis",
                          ncol = NA,
                          size = .1,
                          w_enlarge = 1,
                          h_enlarge = 1
                        ){

                   meta_data   = Seurat@meta.data %>% dplyr::select(target)
                   if(is.na(target_list)%>%all()){target_list = meta_data[,1] %>% levels() %>% as.character()}
                   meta_data_2 = matrix(NA,ncol=length(target_list),nrow=nrow(meta_data)) %>% as.data.frame()
                   rownames(meta_data_2) = rownames(meta_data)
                   colnames(meta_data_2) = target_list
                   
                   for(tmp_target in target_list){
                     meta_data_2[,tmp_target] = ifelse( (meta_data[,1]%>%as.character())==tmp_target,1,0)
                   }

                   colnames(meta_data_2) = colnames(meta_data_2) %>% gsub("-","_",.)
                   target_list = target_list %>% gsub("-","_",.)
                   
                   Seurat_tmp = AddMetaData(Seurat,meta_data_2)
                   
                   p=list()
                   for(iii in 1:length(target_list)){
                       tmp_target   = target_list[iii]
                       p[[iii]] = Nebulosa::plot_density(Seurat_tmp, tmp_target,
                                              reduction = UMAP,
                                              size=size,shape=16,pal=palette) +
                                   theme(legend.position="none")
                       print(tmp_target)
                   }

                   if(is.na(ncol)){ncol_tmp = sqrt(length(target_list) * 2 ) %>% ceiling(.)}else{ncol_tmp=ncol}
                   nrow_tmp = (length(target_list)/ncol_tmp) %>% ceiling(.)
                   
                   p_sum = patchwork::wrap_plots(p, ncol = ncol_tmp) +
                            plot_annotation(paste0(SAMPLE_TMP," : ",name))
                   
                   if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow_tmp*h_enlarge,w=5*ncol_tmp*w_enlarge)}
                   if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow_tmp*h_enlarge,w=480*ncol_tmp*w_enlarge)}
                    plot(p_sum)
                   dev.off()
                   }

######################################################
################# get denstiy score ###############
######################################################

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

get_density_2 <- function(x, y, n=100) {
  h_x  <- bandwidth.nrd(x)
  h_y  <- bandwidth.nrd(y)
  
  if((h_x<=0)){h_x_2 = 0.01}else{h_x_2 = h_x}
  if((h_y<=0)){h_y_2 = 0.01}else{h_y_2 = h_y}
  if( (h_x>0)&(h_y>0) ){
       dens <- MASS::kde2d(x, y, n=n, lims = c(range(x), range(y)))
    }else{
       dens <- MASS::kde2d(x, y, h=c(h_x_2,h_y_2), n=n, lims = c(range(x), range(y)) ) 
    }

  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

######################################################
################# plot split           ###############
######################################################
plot_split = function(Seurat = data,
                        UMAP   = "wnn.RNA.ADT.umap",
                        ident  = "HTO_final_2",
                        target_list = HTO_list_2$HTO,
                        color_list  = NA,
                        with_density = TRUE, 
                        OUTPUT = OUTPUT,
                        name    = name,
                        pt.size = .2,
                        ncol = NA,
                        w_enlarge = 1,
                        h_enlarge = 1
                        ){
                   suppressMessages(library(RColorBrewer))
                   suppressMessages(library(tidyverse))
                   suppressMessages(library(ggrastr))

                   Idents(Seurat)      = ident
                   
                   UMAPdata = Embeddings(Seurat, reduction = UMAP) %>% as.data.frame() %>% rownames_to_column("CellBarcode") %>%
                               left_join(.,Seurat@meta.data %>% as.data.frame() %>% rownames_to_column("CellBarcode"),by="CellBarcode")
                   colnames(UMAPdata)[2:3] = c("UMAP_1","UMAP_2")
                   
                   x_min = min(UMAPdata$UMAP_1) * 1.1
                   x_max = max(UMAPdata$UMAP_1) * 1.1
                   y_min = min(UMAPdata$UMAP_2) * 1.1
                   y_max = max(UMAPdata$UMAP_2) * 1.1
                   
                   ##########################
                   
                   p=list()
                   for(iii in 1:length(target_list)){
                       tmp_target   = target_list[iii]
                       if(is.na(color_list[1])){color="red"}else{color=color_list[iii]}
                       p[[iii]] = ggplot(UMAPdata, aes(x=UMAP_1,y=UMAP_2))
                       if(with_density){p[[iii]] = p[[iii]] + stat_density_2d(geom = "polygon",alpha=0.3,fill="#FEE8C8",color="gray60",breaks=0.002)}
                       p[[iii]] = p[[iii]] +
                                   geom_point_rast(size=.1,color="gray")+
                                   geom_point_rast(data =UMAPdata[UMAPdata[,ident]==tmp_target,],
                                                   aes(x=UMAP_1,y=UMAP_2), size=pt.size,color=color)+
                                   theme_classic()+
                                   xlim(c(x_min,x_max))+
                                   ylim(c(y_min,y_max))+
                                   theme(legend.position="none")+
                                   ggtitle(tmp_target)
                   }
                   
                   if(is.na(ncol)){ncol_tmp = sqrt(length(target_list) * 2 ) %>% ceiling(.)}else{ncol_tmp=ncol}
                   nrow_tmp = (length(target_list)/ncol_tmp) %>% ceiling(.)
                   
                   p_sum = patchwork::wrap_plots(p, ncol = ncol_tmp)+
                            plot_annotation(paste0(SAMPLE_TMP," : ",name," : ","HTO"))
                   
                   ## p[[1]] = AugmentPlot(plot=p[[1]],width=4,height=4,dpi=500)
                   
                   if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5*nrow_tmp*h_enlarge,w=5*ncol_tmp*w_enlarge)}
                   if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=480*nrow_tmp*h_enlarge,w=480*ncol_tmp*w_enlarge)}
                    plot(p_sum)
                   dev.off()
                   }

######################################################
################# plot prop_barplot    ###############
######################################################
plot_prop_barplot = function(Seurat = CITE,
                             x.axis = "tmp",
                             y.axis = "HTO_final_2",
                             theme  = "tmp",
                             OUTPUT = "tmp_porp_barPlot.pdf",
                             col = NA,
                             x_order = NA,
                             y_order = NA,
                             legend = FALSE,
                             col_enlarge = 3,
                             row_enlarge = 1,
                             prop = TRUE
                             ){
                               res_1 = table(Seurat@meta.data[,x.axis],Seurat@meta.data[,y.axis]) %>% as.data.frame.matrix(.)
                               res_2 = res_1 %>% rownames_to_column("x.axis") %>% 
                                        pivot_longer(, col = -x.axis, names_to = "y.axis", values_to = "value") %>% 
                                        as.data.frame()
                               if(is.na(x_order[1])){
                                res_2$x.axis = factor(res_2$x.axis, levels = rownames(res_1))
                               }else{
                                res_2$x.axis = factor(res_2$x.axis, levels = x_order)
                               }

                               if(is.na(y_order[1])){
                                res_2$y.axis = factor(res_2$y.axis, levels = colnames(res_1))
                               }else{
                                res_2$y.axis = factor(res_2$y.axis, levels = y_order)
                               }
                               
                               if(prop==TRUE){p1 = ggplot(res_2,aes(x = x.axis, y = value))+
                                                   geom_bar(stat = "identity", fill = "black" ,color="black")+
                                                   labs(title=theme,y="#cell") +
                                                   theme_void()+
                                                   theme(axis.text.x  = element_blank(),
                                                         axis.text.y  = element_text(colour = "black", size = 10, face = "bold", family = "Helvetica"),
                                                         axis.title.x = element_blank(),
                                                         axis.title.y = element_text(colour = "black", size = 10, face = "bold", family = "Helvetica", angle = 90),
                                                         plot.title   = element_text(colour = "black", size = 10, face = "bold", family = "Helvetica"),
                                                         legend.position="none")
                                             
                                             p2 = ggplot(res_2,aes(x = x.axis, y = value, fill = y.axis))+
                                                   geom_bar(stat = "identity", position="fill",color="black")+
                                                   scale_y_continuous(labels = scales::percent)+
                                                   labs(x="",y="") +
                                                   plot_theme(size=10,theme="void",legend=legend,x_angle=90)
                                             if(!is.na(col[1])){p2 = p2+scale_fill_manual(values=col)}
                                             if(!is.na(x_order[1])){p1 = p1 + scale_x_discrete(limits=x_order)
                                                                    p2 = p2 + scale_x_discrete(limits=x_order)}
              
                                             p_sum = p1 + p2 + plot_layout(ncol = 1, nrow=2, heights = c(.1,.9))
                                }else{
                                             p_sum = ggplot(res_2,aes(x = x.axis, y = value, fill = y.axis))+
                                                      geom_bar(stat = "identity",color="black")+
                                                      plot_theme(size=10,theme="void",legend=legend,x_angle=90)+
                                                      labs(title=theme,x="",y="")
                                             if(!is.na(col[1])){p_sum = p_sum+scale_fill_manual(values=col)}
                                }

                               plot_2(OUTPUT,nrow=row_enlarge,ncol=col_enlarge)
                                plot(p_sum)
                               dev.off()
                              }

######################################################
################# plot differential volcano plot  ####
######################################################

plot_volcano_plot = function(name_tmp   = "name_tmp",
                             Path_tmp   = "differentail_res.txt",
                             OUTPUT     = "Differential_volcano.pdf",
                             log2_FC    = 0.5,
                             pValue     = 0.05,
                             num_repel  = 30,
                             size       = 10,
                             is.repel   = TRUE,
                             repel.size = 1
                             ){ df = fread_FT(Path_tmp)
                                colnames(df)[c(1,3,6)] = c("label","avg_log2FC","p_val_adj")
                              
                                df$log10pVal = -log10(df$p_val_adj)
                                df$log10pVal[df$p_val_adj==0] = 150
                                df$log10pVal[df$log10pVal>150] = 150
                              
                                df$label_2 = NA
                                df$label_2[1:num_repel] = df$label[1:num_repel]
                              
                                df$DEG = "nonDEG"
                                df$DEG[(df$avg_log2FC>log2_FC)&(df$p_val_adj<pValue)] = "posiDEG"
                                df$DEG[(-log2_FC>df$avg_log2FC)&(df$p_val_adj<pValue)] = "negaDEG"
                                df$DEG = df$DEG %>% factor(.,levels=c("posiDEG","negaDEG","nonDEG"))

                                posi_DEG_num = sum(df$DEG == "posiDEG")
                                nega_DEG_num = sum(df$DEG == "negaDEG")
                              
                                col_DEG = data.frame(DEG=c("posiDEG","negaDEG","nonDEG"),
                                                   color=c("red","blue","white"))
                                col_DEG_2 = col_DEG[unique(df$DEG) %>% as.character() %>% is.element(col_DEG$DEG,.),]
                              
                                p1 = ggplot(df,aes(x=avg_log2FC,y=log10pVal,fill=DEG,label=label_2))+
                                      geom_vline(xintercept=log2_FC,col="grey",size=1,linetype="dashed")+
                                      geom_vline(xintercept=-log2_FC,col="grey",size=1,linetype="dashed")+
                                      geom_hline(yintercept=-log10(pValue),col="grey",size=1,linetype="dashed")+
                                      geom_point(size=1,shape=21,color="black")+
                                      plot_theme(size=size,legend=FALSE)+
                                      scale_fill_manual(values=col_DEG_2$color)+
                                      labs(title=paste0(name_tmp," : differential : posi ",posi_DEG_num," : nega ",nega_DEG_num),
                                           x="log2 FC",y="-log10(adjusted_pValue)",fill="DEG")
                                      
                                if(is.repel){p1 = p1 + ggrepel::geom_text_repel(col="black",fontface = "bold",size=repel.size,na.rm=TRUE) }
                              
                                if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5,w=5)}
                                if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=960,w=960)}
                                 plot(p1)
                                dev.off()
}

plot_volcano_rast = function(df = peak_diff,
                             name_tmp,
                             perc_filter= 0.05,
                             log2_FC    = 0.25,
                             log2_FC_thresh = 3,
                             pValue     = 0.05,
                             num_repel  = 30,
                             size       = 10,
                             is.repel   = TRUE,
                             repel.size = 1
                        ){ colnames(df)[c(1,4,5)] = c("label","perc_trgt","perc_ctrl")
                           df = df %>%
                                 dplyr::filter((perc_trgt>=perc_filter)|(perc_ctrl>=perc_filter)) %>%
                                 mutate(avg_log2FC = case_when((avg_log2FC>=log2_FC_thresh)~log2_FC_thresh,TRUE~avg_log2FC))
                           
                           df$log10pVal = -log10(df$p_val_adj)
                           df$log10pVal[df$p_val_adj==0] = 300
                           df$log10pVal[df$log10pVal>300] = 300
                         
                           df$label_2 = NA
                           df$label_2[1:num_repel] = df$label[1:num_repel]
                         
                           df$DEG = "nonDEG"
                           df$DEG[(df$avg_log2FC>log2_FC)&(df$p_val_adj<pValue)] = "posiDEG"
                           df$DEG[(-log2_FC>df$avg_log2FC)&(df$p_val_adj<pValue)] = "negaDEG"
                           df$DEG = df$DEG %>% factor(.,levels=c("posiDEG","negaDEG","nonDEG"))

                           posi_DEG_num = sum(df$DEG == "posiDEG")
                           nega_DEG_num = sum(df$DEG == "negaDEG")
                         
                           col_DEG = data.frame(DEG=c("posiDEG","negaDEG","nonDEG"),
                                              color=c("red","blue","white"))
                           col_DEG_2 = col_DEG[unique(df$DEG) %>% as.character() %>% is.element(col_DEG$DEG,.),]
                         
                           p1 = ggplot(df,aes(x=avg_log2FC,y=log10pVal,fill=DEG,label=label_2))+
                                 geom_vline(xintercept=log2_FC,col="grey",linewidth=1,linetype="dashed")+
                                 geom_vline(xintercept=-log2_FC,col="grey",linewidth=1,linetype="dashed")+
                                 geom_hline(yintercept=-log10(pValue),col="grey",linewidth=1,linetype="dashed")+
                                 geom_point_rast(size=1,shape=21,color="black")+
                                 plot_theme(theme="void",size=size,legend=FALSE)+
                                 scale_fill_manual(values=col_DEG_2$color)+
                                 scale_x_continuous(breaks=c(-1*log2_FC_thresh,0,log2_FC_thresh),limits=c(-1.1*log2_FC_thresh,1.1*log2_FC_thresh))+
                                 scale_y_continuous(breaks=c(0,150,300),limits=c(-10,300))+
                                 labs(title=paste0(name_tmp," : differential : posi ",posi_DEG_num," : nega ",nega_DEG_num),
                                      x="log2 FC",y="-log10(adjusted_pValue)",fill="DEG")
                                 
                           if(is.repel){p1 = p1 + ggrepel::geom_text_repel(col="black",fontface = "bold",size=repel.size,na.rm=TRUE) }
                           p1
}

######################################################
################# plot differential MA plot  ####
######################################################

plot_MA_plot = function(name_tmp   = "name_tmp",
                        Path_tmp   = "differentail_res.txt",
                        OUTPUT     = "Differential_volcano.pdf",
                        num_repel  = 30,
                        size       = 10,
                        is.repel   = TRUE,
                        repel.size = 1
                             ){ df = fread_FT(Path_tmp)
                                df$DEG = "nonDEG"
                                df$DEG[(df$FDR<0.05)&(df$logFC>0)] = "posiDEG"
                                df$DEG[(df$FDR<0.01)&(df$logFC>0)] = "posiDEG_2"
                                df$DEG[(df$FDR<0.05)&(df$logFC<0)] = "negaDEG"
                                df$DEG[(df$FDR<0.01)&(df$logFC<0)] = "negaDEG_2"
                                df$label = NA
                                df$label[1:num_repel] = df$Gene[1:num_repel]
                                
                                posi_DEG_num   = sum(df$DEG == "posiDEG_2")
                                nega_DEG_num   = sum(df$DEG == "negaDEG_2")
                                posi_DEG_num_2 = sum(df$DEG %in% c("posiDEG","posiDEG_2") )
                                nega_DEG_num_2 = sum(df$DEG %in% c("negaDEG","negaDEG_2") )
                                
                                col_DEG = data.frame(DEG=c("posiDEG_2","posiDEG","negaDEG_2","negaDEG","nonDEG"),
                                                     color=c("red","#FF7F7F","blue","#7F7FFF","white"))
                                col_DEG_2 = col_DEG[unique(df$DEG) %>% as.character() %>% is.element(col_DEG$DEG,.),]
                                df$DEG = factor(df$DEG,levels=col_DEG_2$DEG)
                                
                                p1 = ggplot(df,aes(x=logCPM,y=logFC,fill=DEG,label=label))+
                                      geom_hline(yintercept=1,col="gray",size=.5)+
                                      geom_hline(yintercept=.5,col="gray",size=.5)+
                                      geom_hline(yintercept=0,col="black",size=1)+
                                      geom_hline(yintercept=-.5,col="gray",size=.5)+
                                      geom_hline(yintercept=-1,col="gray",size=.5)+
                                      geom_point(size=1,shape=21,color="black")+
                                      plot_theme(size=size,legend=FALSE)+
                                      scale_fill_manual(values=col_DEG_2$color)+
                                      labs(title=paste0(name_tmp," : differential : posi FDR0.05 ",posi_DEG_num_2," FDR0.01 ",posi_DEG_num," : nega FDR0.05 ",nega_DEG_num_2," FDR0.01 ",nega_DEG_num),
                                           x="Average logCPM",y="logFC" )
                                      
                                if(is.repel){p1 = p1 + ggrepel::geom_text_repel(col="black",fontface = "bold",size=repel.size,na.rm=TRUE) }

                                if(grepl("pdf$",OUTPUT)){pdf_2(OUTPUT,h=5,w=5)}
                                if(grepl("png$",OUTPUT)){png_2(OUTPUT,h=960,w=960)}
                                 plot(p1)
                                dev.off()
}

######################################################
################# Assay quantile filter     ##########
######################################################

Assay_quant_file = function(Seurat = Seurat,
                            assay  = "ADT",
                            assay_filt = "ADTQuan",
                            quantile = 0.99){
                                             tmp_counts = GetAssayData(Seurat, assay = assay , slot = "counts")
                                             Quan_assay = t(tmp_counts) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,quantile)})
                                             omit_tmp   = names(Quan_assay)[Quan_assay==0]
                                             Seurat[[assay_filt]] = CreateAssayObject(counts = tmp_counts[!is.element(rownames(tmp_counts),omit_tmp),] )
                                             Seurat
                                            }

######################################################
################# Assay quantile filter     ########## 230615訂正済
######################################################

Assay_quant_data = function(Seurat = Seurat,
                            assay  = "ADT",
                            assay_filt = "ADTquan",
                            quantile = 0.99,
                            low_qunatile = 0,
                            sep=1          ){if(sep==1){####
                                                        tmp_counts = GetAssayData(Seurat, assay = assay , slot = "counts")
                                                        Quan_assay = t(tmp_counts) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,quantile)})
                                                        Quan_assay_l = t(tmp_counts) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,low_quantile)})
                                                        omit_tmp   = names(Quan_assay)[Quan_assay==Quan_assay_l]
                                                        tmp_counts = tmp_counts[!is.element(rownames(tmp_counts),omit_tmp),]
                                                        Quan_assay = Quan_assay[rownames(tmp_counts)]
                                                        Quan_assay_l = Quan_assay_l[rownames(tmp_counts)]
                                                        for(tmp_index in 1:length(Quan_assay)){
                                                          tmp_counts[tmp_index,] = tmp_counts[tmp_index,] %>%
                                                                                    sapply(.,function(x){x[x>=Quan_assay[tmp_index]]=Quan_assay[tmp_index];x}) %>%
                                                                                    sapply(.,function(x){x[x<=Quan_assay_l[tmp_index]]=Quan_assay_l[tmp_index];x})
                                                          show.progress(tmp_index,1:length(Quan_assay))
                                                        }
                                                        Seurat[[assay_filt]] = CreateAssayObject(counts = tmp_counts )
                                                        ####
                                                        tmp_data = GetAssayData(Seurat, assay = assay , slot = "data")
                                                        Quan_assay = t(tmp_data) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,quantile)})
                                                        Quan_assay_l = t(tmp_data) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,low_quantile)})
                                                        omit_tmp   = names(Quan_assay)[Quan_assay==Quan_assay_l]
                                                        tmp_data = tmp_data[!is.element(rownames(tmp_data),omit_tmp),]
                                                        Quan_assay = Quan_assay[rownames(tmp_data)]
                                                        Quan_assay_l = Quan_assay_l[rownames(tmp_data)]
                                                        for(tmp_index in 1:length(Quan_assay)){
                                                          tmp_data[tmp_index,] = tmp_data[tmp_index,] %>%
                                                                                  sapply(.,function(x){x[x>=Quan_assay[tmp_index]]=Quan_assay[tmp_index];x}) %>%
                                                                                  sapply(.,function(x){x[x<=Quan_assay_l[tmp_index]]=Quan_assay_l[tmp_index];x})
                                                          show.progress(tmp_index,1:length(Quan_assay))
                                                        }
                                                        Seurat = SetAssayData(object=Seurat, assay=assay_filt, slot="data", new.data = tmp_data )
                                                        Seurat = ScaleData(Seurat, assay=assay_filt)
                                                       }
                                             if(sep>1 ){####
                                                        tmp_counts = GetAssayData(Seurat, assay = assay , slot = "counts")
                                                        each_num = ceiling( nrow(tmp_counts) / sep )
                                                        tmp_counts_list = list()
                                                        for(tmp_num in 1:sep){tmp_from = each_num * (tmp_num - 1) + 1
                                                                              tmp_to   = min(each_num * tmp_num , nrow(tmp_counts) )
                                                                              tmp_counts_tmp = tmp_counts[tmp_from:tmp_to,]
                                                                              Quan_assay = t(tmp_counts_tmp) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,quantile)})
                                                                              Quan_assay_l = t(tmp_counts_tmp) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,low_quantile)})
                                                                              omit_tmp   = names(Quan_assay)[Quan_assay==Quan_assay_l]
                                                                              tmp_counts_tmp = tmp_counts_tmp[!is.element(rownames(tmp_counts_tmp),omit_tmp),]
                                                                              Quan_assay = Quan_assay[rownames(tmp_counts_tmp)]
                                                                              Quan_assay_l = Quan_assay_l[rownames(tmp_counts_tmp)]
                                                                              for(tmp_index in 1:length(Quan_assay)){
                                                                                  tmp_counts_tmp[tmp_index,] = tmp_counts_tmp[tmp_index,] %>%
                                                                                                                sapply(.,function(x){x[x>=Quan_assay[tmp_index]]=Quan_assay[tmp_index];x}) %>%
                                                                                                                sapply(.,function(x){x[x<=Quan_assay_l[tmp_index]]=Quan_assay_l[tmp_index];x})
                                                                                  ## show.progress(tmp_index,1:length(Quan_assay))
                                                                                 }
                                                                              tmp_counts_list[[tmp_num]] = tmp_counts_tmp
                                                                              show.progress(tmp_num,1:sep)
                                                                              }
                                                        Seurat[[assay_filt]] = CreateAssayObject(counts = do.call(rbind,tmp_counts_list) )
                                                        ####
                                                        tmp_data = GetAssayData(Seurat, assay = assay , slot = "data")
                                                        each_num = ceiling( nrow(tmp_data) / sep )
                                                        tmp_data_list = list()
                                                        for(tmp_num in 1:sep){tmp_from = each_num * (tmp_num - 1) + 1
                                                                              tmp_to   = min(each_num * tmp_num , nrow(tmp_data) )
                                                                              tmp_data_tmp = tmp_data[tmp_from:tmp_to,]
                                                                              Quan_assay = t(tmp_data_tmp) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,quantile)})
                                                                              Quan_assay_l = t(tmp_data_tmp) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,low_quantile)})
                                                                              omit_tmp   = names(Quan_assay)[Quan_assay==Quan_assay_l]
                                                                              tmp_data_tmp = tmp_data_tmp[!is.element(rownames(tmp_data_tmp),omit_tmp),]
                                                                              Quan_assay = Quan_assay[rownames(tmp_data_tmp)]
                                                                              Quan_assay_l = Quan_assay_l[rownames(tmp_counts_tmp)]
                                                                              for(tmp_index in 1:length(Quan_assay)){
                                                                                  tmp_data_tmp[tmp_index,] = tmp_data_tmp[tmp_index,] %>%
                                                                                                              sapply(.,function(x){x[x>=Quan_assay[tmp_index]]=Quan_assay[tmp_index];x}) %>%
                                                                                                              sapply(.,function(x){x[x<=Quan_assay_l[tmp_index]]=Quan_assay_l[tmp_index];x})
                                                                                  ## show.progress(tmp_index,1:length(Quan_assay))
                                                                                 }
                                                                              tmp_data_list[[tmp_num]] = tmp_data_tmp
                                                                              show.progress(tmp_num,1:sep)
                                                                              }
                                                        Seurat = SetAssayData(object=Seurat, assay=assay_filt, slot="data", new.data = do.call(rbind,tmp_data_list) )
                                                        Seurat = ScaleData(Seurat, assay=assay_filt)
                                                       }

                                             ####
                                             Seurat
                                            }


######################################################
################# CITE :: SCT normalize     ##########
######################################################

CITE_SCT_WNN = function(CITE = CITE,
                        assay_RNA       = "RNA",
                        var_to_regress  = c("percent.mt"),
                        with.RNA.filt   = FALSE,
                        assay_RNA_filt  = "RNAfilt",
                        min_RNA_filt    = 0, 
                        assay_RNA_SCT   = "RNA.SCT",
                        residual.features = NULL, 
                        red_RNA_SCT     = "RNA.SCT.pca",
                        redkey_RNA_SCT  = "RNASCTpca_",
                        graph_RNA_SCT   = "RNA.SCT.pca.nn",
                        red2_RNA_SCT    = "RNA.SCT.pca.umap",
                        redkey2_RNA_SCT = "RNASCTpcaUMAP_",
                        dim_RNA_SCT     = 1:30,
                        res_RNA_SCT     = 0.8,
                        assay_ADT       = "ADT",
                        assay_ADT_filt  = "ADTfilt",
                        ADT_filt_thresh = 0.75,
                        red_ADT         = "ADT.pca",
                        redkey_ADT      = "ADTpca_",
                        graph_ADT       = "ADT.pca.nn",
                        red2_ADT        = "ADT.pca.umap",
                        redkey2_ADT     = "ADTpcaUMAP_",
                        dim_ADT         = 2:30,
                        res_ADT         = 0.8,
                        WNN.modality.weight.name = c("RNA.weight","ADT.weight"),
                        WNN.weighted.nn.name = "weighted.RNA.ADT.nn",
                        WNN.knn.graph.name = "RNA.ADT.wknn",
                        WNN.snn.graph.name = "RNA.ADT.wsnn",
                        red_WNN         = "wnn.RNA.ADT.umap",
                        redkey_WNN      = "wnnRNAADTUMAP_",
                        algorithm_WNN   = 3,
                        res_WNN         = 0.8,
                        omit_gene_PATH  = NA, gene_for_use = NA
                        ){
                         #####################################################################
                         #### scRNA UMAP
                         #####################################################################
                         if(with.RNA.filt){RNA_counts = GetAssayData(CITE, assay = assay_RNA , slot = "counts") 
                                           CITE[[assay_RNA_filt]] = CreateAssayObject(counts = RNA_counts[rowSums(RNA_counts)>min_RNA_filt,] )
                                           CITE = SCTransform(CITE,assay = assay_RNA_filt , new.assay.name = assay_RNA_SCT, vars.to.regress = var_to_regress, residual.features = residual.features, verbose = FALSE, return.only.var.genes = FALSE)}
                         if(!with.RNA.filt){CITE = SCTransform(CITE,assay = assay_RNA , new.assay.name = assay_RNA_SCT, vars.to.regress = var_to_regress, residual.features = residual.features, verbose = FALSE, return.only.var.genes = FALSE)}
                         ## CITE = ScaleData(CITE , assay = assay_RNA_SCT) ## assay = assay_RNA
                         variable_gene = VariableFeatures(object = CITE, assay = assay_RNA_SCT )
                         if(!is.na(omit_gene_PATH)){omit_gene = read.table_FT(omit_gene_PATH)[,1];variable_gene = setdiff(variable_gene,omit_gene)}
                         if(!is.na(gene_for_use[1])){variable_gene = gene_for_use}
                         CITE = RunPCA(CITE, assay = assay_RNA_SCT , reduction.name = red_RNA_SCT , reduction.key = redkey_RNA_SCT, verbose = FALSE,
                                       features = variable_gene )
                         CITE = FindNeighbors(CITE, assay=assay_RNA_SCT, reduction = red_RNA_SCT , graph.name = graph_RNA_SCT, dims = dim_RNA_SCT)
                         CITE = FindClusters(CITE, graph.name = graph_RNA_SCT , resolution = res_RNA_SCT )
                         CITE = RunUMAP(CITE, dims = dim_RNA_SCT, reduction = red_RNA_SCT ,reduction.name = red2_RNA_SCT, reduction.key = redkey2_RNA_SCT)
                         
                         #####################################################################
                         #### protein data UMAP
                         #####################################################################
                         
                         ADT_counts = GetAssayData(CITE, assay = assay_ADT , slot = "counts")
                         ADT_data   = GetAssayData(CITE, assay = assay_ADT , slot = "data")
                         Quan_ADT = t(ADT_counts) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,ADT_filt_thresh)})
                         omit_ADT      = names(Quan_ADT)[Quan_ADT==0]
                         
                         CITE[[assay_ADT_filt]] = CreateAssayObject(counts = ADT_counts[!is.element(rownames(ADT_counts),omit_ADT),] )
                         CITE                   = SetAssayData(object=CITE, assay = assay_ADT_filt , slot="data", 
                                                               new.data = ADT_data[!is.element(rownames(ADT_data),omit_ADT),] %>% as.matrix(.) )
                         CITE                   = ScaleData(CITE , assay = assay_ADT_filt)
                         
                         DefaultAssay(CITE)     = assay_ADT_filt
                         variable_feature = rownames(CITE[[assay_ADT_filt]])
                         CITE = RunPCA(CITE   , assay=assay_ADT_filt, reduction.name = red_ADT, reduction.key = redkey_ADT , verbose = FALSE, features = variable_feature )
                         CITE = FindNeighbors(CITE, assay = assay_ADT_filt, reduction = red_ADT , graph.name = graph_ADT , dims = dim_ADT)
                         CITE = FindClusters(CITE, graph.name = graph_ADT , resolution = res_ADT)
                         CITE = RunUMAP(CITE, dims = dim_ADT , reduction = red_ADT , reduction.name = red2_ADT , reduction.key = redkey2_ADT)
                         
                         #####################################################################
                         #### WNN -- scRNA and ADT
                         #####################################################################
                         
                         CITE = FindMultiModalNeighbors(CITE, 
                                                              reduction.list       = list(red_RNA_SCT, red_ADT), 
                                                              dims.list            = list(dim_RNA_SCT, dim_ADT), 
                                                              modality.weight.name = WNN.modality.weight.name,
                                                              weighted.nn.name     = WNN.weighted.nn.name,
                                                              knn.graph.name       = WNN.knn.graph.name,
                                                              snn.graph.name       = WNN.snn.graph.name
                                                              )
                         CITE = RunUMAP(CITE, nn.name = WNN.weighted.nn.name, reduction.name = red_WNN, reduction.key = redkey_WNN)
                         CITE = FindClusters(CITE, graph.name = WNN.snn.graph.name , algorithm = algorithm_WNN, resolution = res_WNN, verbose = FALSE)

                         CITE
                    }

######################################################
################# CITE :: SCT normalize w/ cellcyle ##
######################################################

CITE_SCT_CellCycle_WNN = function(CITE = CITE,
                                  assay_RNA       = "RNA",
                                  var_to_regress  = c("percent.mt","S.Score","G2M.Score"),
                                  with.RNA.filt   = FALSE,
                                  assay_RNA_filt  = "RNAfilt",
                                  min_RNA_filt    = 0, 
                                  assay_RNA_SCT   = "RNA.SCT.CC",
                                  residual.features = NULL,
                                  red_RNA_SCT     = "RNA.SCT.CC.pca",
                                  redkey_RNA_SCT  = "RNASCTCCpca_",
                                  graph_RNA_SCT   = "RNA.SCT.CC.pca.nn",
                                  red2_RNA_SCT    = "RNA.SCT.CC.pca.umap",
                                  redkey2_RNA_SCT = "RNASCTCCpcaUMAP_",
                                  dim_RNA_SCT     = 1:30,
                                  res_RNA_SCT     = 0.8,
                                  assay_ADT       = "ADT",
                                  assay_ADT_filt  = "ADTfilt",
                                  ADT_filt_thresh = 0.75,
                                  red_ADT         = "ADT.pca",
                                  redkey_ADT      = "ADTpca_",
                                  graph_ADT       = "ADT.pca.nn",
                                  red2_ADT        = "ADT.pca.umap",
                                  redkey2_ADT     = "ADTpcaUMAP_",
                                  dim_ADT         = 2:30,
                                  res_ADT         = 0.8,
                                  WNN.modality.weight.name = c("RNA.weight","ADT.weight"),
                                  WNN.weighted.nn.name = "weighted.RNA.ADT.CC.nn",
                                  WNN.knn.graph.name = "RNA.ADT.CC.wknn",
                                  WNN.snn.graph.name = "RNA.ADT.CC.wsnn",
                                  red_WNN         = "wnn.RNA.ADT.CC.umap",
                                  redkey_WNN      = "wnnRNAADTCCUMAP_",
                                  algorithm_WNN   = 3,
                                  res_WNN         = 0.8,
                                  omit_gene_PATH  = NA, gene_for_use = NA
                                  ){
                                   #####################################################################
                                   #### scRNA UMAP
                                   #####################################################################
                                   if(with.RNA.filt){RNA_counts = GetAssayData(CITE, assay = assay_RNA , slot = "counts") 
                                                     CITE[[assay_RNA_filt]] = CreateAssayObject(counts = RNA_counts[rowSums(RNA_counts)>min_RNA_filt,] )
                                                     CITE = SCTransform(CITE,assay = assay_RNA_filt , new.assay.name = assay_RNA_SCT, vars.to.regress = var_to_regress, residual.features = residual.features, verbose = FALSE, return.only.var.genes = FALSE)}
                                   if(!with.RNA.filt){CITE = SCTransform(CITE,assay = assay_RNA , new.assay.name = assay_RNA_SCT, vars.to.regress = var_to_regress, residual.features = residual.features, verbose = FALSE, return.only.var.genes = FALSE)}
                                   ## CITE = ScaleData(CITE , assay = assay_RNA)
                                   variable_gene = VariableFeatures(object = CITE, assay = assay_RNA_SCT )
                                   if(!is.na(omit_gene_PATH)){omit_gene = read.table_FT(omit_gene_PATH)[,1];variable_gene = setdiff(variable_gene,omit_gene)}
                                   if(!is.na(gene_for_use[1])){variable_gene = gene_for_use}
                                   CITE = RunPCA(CITE, assay = assay_RNA_SCT , reduction.name = red_RNA_SCT , reduction.key = redkey_RNA_SCT, verbose = FALSE,
                                                 features = variable_gene )
                                   CITE = FindNeighbors(CITE, assay=assay_RNA_SCT, reduction = red_RNA_SCT , graph.name = graph_RNA_SCT, dims = dim_RNA_SCT)
                                   CITE = FindClusters(CITE, graph.name = graph_RNA_SCT , resolution = res_RNA_SCT )
                                   CITE = RunUMAP(CITE, dims = dim_RNA_SCT, reduction = red_RNA_SCT ,reduction.name = red2_RNA_SCT, reduction.key = redkey2_RNA_SCT)
                                   
                                   #####################################################################
                                   #### protein data UMAP
                                   #####################################################################
                                   
                                   ADT_counts = GetAssayData(CITE, assay = assay_ADT , slot = "counts")
                                   ADT_data   = GetAssayData(CITE, assay = assay_ADT , slot = "data")
                                   Quan_ADT = t(ADT_counts) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,ADT_filt_thresh)})
                                   omit_ADT      = names(Quan_ADT)[Quan_ADT==0]
                                   
                                   CITE[[assay_ADT_filt]] = CreateAssayObject(counts = ADT_counts[!is.element(rownames(ADT_counts),omit_ADT),] )
                                   CITE                   = SetAssayData(object=CITE, assay = assay_ADT_filt , slot="data", 
                                                                         new.data = ADT_data[!is.element(rownames(ADT_data),omit_ADT),] %>% as.matrix(.) )
                                   CITE                   = ScaleData(CITE , assay = assay_ADT_filt)
                                   
                                   DefaultAssay(CITE)     = assay_ADT_filt
                                   variable_feature = rownames(CITE[[assay_ADT_filt]])
                                   CITE = RunPCA(CITE   , assay=assay_ADT_filt, reduction.name = red_ADT, reduction.key = redkey_ADT , verbose = FALSE, features = variable_feature )
                                   CITE = FindNeighbors(CITE, assay = assay_ADT_filt, reduction = red_ADT , graph.name = graph_ADT , dims = dim_ADT)
                                   CITE = FindClusters(CITE, graph.name = graph_ADT , resolution = res_ADT)
                                   CITE = RunUMAP(CITE, dims = dim_ADT , reduction = red_ADT , reduction.name = red2_ADT , reduction.key = redkey2_ADT)
                                   
                                   #####################################################################
                                   #### WNN -- scRNA and ADT
                                   #####################################################################
                                   
                                   CITE <- FindMultiModalNeighbors(CITE, 
                                                                        reduction.list       = list(red_RNA_SCT, red_ADT), 
                                                                        dims.list            = list(dim_RNA_SCT, dim_ADT), 
                                                                        modality.weight.name = WNN.modality.weight.name,
                                                                        weighted.nn.name     = WNN.weighted.nn.name,
                                                                        knn.graph.name       = WNN.knn.graph.name,
                                                                        snn.graph.name       = WNN.snn.graph.name
                                                                        )
                                   CITE = RunUMAP(CITE, nn.name = WNN.weighted.nn.name, reduction.name = red_WNN, reduction.key = redkey_WNN)
                                   CITE = FindClusters(CITE, graph.name = WNN.snn.graph.name , algorithm = algorithm_WNN, resolution = res_WNN, verbose = FALSE)

                                   CITE
                                }


######################################################
################# ASAP normalize            ##########
######################################################

ASAP_WNN = function(ASAP = ASAP,
                    assay_ATAC      = "ATAC",
                    red_ATAC        = "ATAC.lsi",
                    redkey_ATAC     = "ATACLSI_",
                    graph_ATAC      = "ATAC.lsi.nn",
                    red2_ATAC       = "ATAC.lsi.umap",
                    redkey2_ATAC    = "ATAClsiUMAP_",
                    dim_ATAC        = 2:30,
                    res_ATAC        = 0.8,
                    algorithm_ATAC  = 3,
                    assay_ADT       = "ADT",
                    assay_ADT_filt  = "ADTfilt",
                    ADT_filt_thresh = 0.75,
                    red_ADT         = "ADT.pca",
                    redkey_ADT      = "ADTpca_",
                    graph_ADT       = "ADT.pca.nn",
                    red2_ADT        = "ADT.pca.umap",
                    redkey2_ADT     = "ADTpcaUMAP_",
                    dim_ADT         = 2:30,
                    res_ADT         = 0.8,
                    WNN.weighted.nn.name = "weighted.ATAC.ADT.nn",
                    WNN.knn.graph.name = "ATAC.ADT.wknn",
                    WNN.snn.graph.name = "ATAC.ADT.wsnn",
                    red_WNN         = "wnn.ATAC.ADT.umap",
                    redkey_WNN      = "wnnATACADTUMAP_",
                    algorithm_WNN   = 3,
                    res_WNN         = 0.8
                        ){
                          #####################################################################
                          #### ATAC peak --- TFIDF UMP 
                          #####################################################################
                          
                          DefaultAssay(ASAP) <- assay_ATAC
                          #Dimensional Reduction
                          ASAP = RunTFIDF(ASAP, assay = assay_ATAC)
                          ASAP = FindTopFeatures(ASAP, assay = assay_ATAC , min.cutoff = 'q0')
                          ASAP = RunSVD(ASAP, assay = assay_ATAC ,reduction.name = red_ATAC , reduction.key = redkey_ATAC )
                          
                          # UMAP
                          ASAP  = RunUMAP(object = ASAP, reduction = red_ATAC , dims = dim_ATAC ,reduction.name = red2_ATAC , reduction.key = redkey2_ATAC )
                          ASAP  = FindNeighbors(object = ASAP, reduction = red_ATAC , graph.name = graph_ATAC , dims = dim_ATAC)
                          ASAP  = FindClusters(object = ASAP,  graph.name = graph_ATAC , verbose = FALSE, algorithm = algorithm_ATAC , resolution = res_ATAC)
                          
                          #####################################################################
                          #### protein data UMAP
                          #####################################################################
                          
                          ADT_counts = GetAssayData(ASAP, assay = assay_ADT , slot = "counts")
                          ADT_data   = GetAssayData(ASAP, assay = assay_ADT , slot = "data")
                          Quan_ADT   = t(ADT_counts) %>% as.data.frame(.) %>% apply(.,2,function(x){quantile(x,ADT_filt_thresh)})
                          omit_ADT   = names(Quan_ADT)[Quan_ADT==0]
                          
                          ASAP[[assay_ADT_filt]] = CreateAssayObject(counts = ADT_counts[!is.element(rownames(ADT_counts),omit_ADT),] )
                          ASAP                   = SetAssayData(object=ASAP, assay = assay_ADT_filt , slot = "data" , 
                                                                new.data = ADT_data[!is.element(rownames(ADT_data),omit_ADT),] %>% as.matrix(.) )
                          ASAP                   = ScaleData( ASAP ,    assay = assay_ADT_filt )
                          
                          DefaultAssay(ASAP)     = assay_ADT_filt
                          VariableFeatures(ASAP) = rownames(ASAP[[assay_ADT_filt]])
                          ASAP = RunPCA(ASAP   , assay = assay_ADT_filt, reduction.name = red_ADT , reduction.key = redkey_ADT, verbose = FALSE, features = ASAP[[assay_ADT_filt]] %>% rownames(.) )
                          ASAP = FindNeighbors(ASAP, assay = assay_ADT_filt, reduction = red_ADT , graph.name = graph_ADT , dims = dim_ADT )
                          ASAP = FindClusters(ASAP, graph.name = graph_ADT , resolution = res_ADT )
                          ASAP = RunUMAP(ASAP, dims = dim_ADT , reduction = red_ADT ,reduction.name = red2_ADT , reduction.key = redkey2_ADT )
                          
                          #####################################################################
                          #### WNN -- scRNA and ADT
                          #####################################################################
                          
                          ASAP = FindMultiModalNeighbors(ASAP, 
                                                         reduction.list   = list(red_ATAC, red_ADT), 
                                                         dims.list        = list(dim_ATAC, dim_ADT), 
                                                         weighted.nn.name = WNN.weighted.nn.name,
                                                         knn.graph.name   = WNN.knn.graph.name,
                                                         snn.graph.name   = WNN.snn.graph.name
                                                         )
                          
                          ASAP = RunUMAP(ASAP, nn.name = WNN.weighted.nn.name, reduction.name = red_WNN, reduction.key = redkey_WNN)
                          ASAP = FindClusters(ASAP, graph.name = WNN.snn.graph.name , algorithm = algorithm_WNN, resolution = res_WNN, verbose = FALSE)
                          
                          ASAP
                        }

######################################################
################# ASAP normalize Harmony    ##########
######################################################

ASAP_WNN_harmony = function(ASAP = ASAP,
                    assay_ATAC      = "ATAC",
                    var_to_regress  = c("percent.mt","S.Score","G2M.Score"),
                    red_ATAC        = "ATAC.lsi",
                    red_ATAC_harm   = "ATAC.Harmony",
                    graph_ATAC      = "ATAC.Harmony.nn",
                    red2_ATAC       = "ATAC.Harmony.umap",
                    redkey2_ATAC    = "ATACHarmonyumap_",
                    dim_ATAC        = 2:30,
                    res_ATAC        = 0.8,
                    algorithm_ATAC  = 3,
                    assay_ADT_filt  = "ADTfilt",
                    red_ADT         = "ADT.pca",
                    red_ADT_harm    = "ADT.Harmony",
                    graph_ADT       = "ADT.Harmony.nn",
                    red2_ADT        = "ADT.Harmony.umap",
                    redkey2_ADT     = "ADTHarmonyumap_",
                    dim_ADT         = 2:30,
                    res_ADT         = 0.8,
                    WNN.weighted.nn.name = "weighted.ATAC.ADT.Harmony.nn",
                    WNN.knn.graph.name = "ATAC.ADT.Harmony.wknn",
                    WNN.snn.graph.name = "ATAC.ADT.Harmony.wsnn",
                    red_WNN         = "wnn.ATAC.ADT.Harmony.umap",
                    redkey_WNN      = "wnnATACADTHarmonyUMAP_",
                    algorithm_WNN   = 3,
                    res_WNN         = 0.8
                        ){
                          #####################################################################
                          #### ATAC peak --- TFIDF UMP --- Harmony
                          #####################################################################
                          
                          ASAP = RunHarmony(ASAP, group.by.vars=var_to_regress, reduction = red_ATAC, assay.use = assay_ATAC, reduction.save = red_ATAC_harm, project.dim = FALSE)
                          ASAP = RunUMAP(object = ASAP, reduction = red_ATAC_harm, dims = dim_ATAC, reduction.name = red2_ATAC, reduction.key = redkey2_ATAC)
                          ASAP = FindNeighbors(object = ASAP, reduction = red_ATAC_harm, graph.name = graph_ATAC, dims = dim_ATAC)
                          ASAP = FindClusters(object = ASAP, graph.name = graph_ATAC, verbose = FALSE, resolution = res_ATAC , algorithm = algorithm_ATAC)

                          #####################################################################
                          #### protein data UMAP --- Harmony
                          #####################################################################
                          
                          ASAP = RunHarmony(ASAP, group.by.vars=var_to_regress, reduction = red_ADT, assay.use = assay_ADT_filt, reduction.save = red_ADT_harm)
                          ASAP = RunUMAP(object = ASAP, reduction = red_ADT_harm, dims = dim_ADT, reduction.name = red2_ADT, reduction.key = redkey2_ADT)
                          ASAP = FindNeighbors(object = ASAP, reduction = red_ADT_harm,graph.name = graph_ADT, dims = dim_ADT)
                          ASAP = FindClusters(object = ASAP, graph.name = graph_ADT, resolution = res_ADT )
                          
                          #####################################################################
                          #### WNN -- scRNA and ADT --- Harmony
                          #####################################################################
                          
                          ASAP = FindMultiModalNeighbors(ASAP, 
                                                          reduction.list = list(red_ATAC_harm, red_ADT_harm), 
                                                          dims.list = list(dim_ATAC, dim_ADT),
                                                          weighted.nn.name = WNN.weighted.nn.name,
                                                          knn.graph.name = WNN.knn.graph.name,
                                                          snn.graph.name = WNN.snn.graph.name
                                                          )
                          
                          ASAP = RunUMAP(ASAP, nn.name = WNN.weighted.nn.name, reduction.name = red_WNN, reduction.key = redkey_WNN)
                          ASAP = FindClusters(ASAP, graph.name = WNN.snn.graph.name, algorithm = algorithm_WNN, resolution = res_WNN, verbose = FALSE)
                          
                          ASAP
                        }

######################################################
################# Seurat HTO UAMP      ###############
######################################################

HTO_UMAP = function(Seurat = data,
                    assay      = "HTO",
                    red        = "HTO.pca",
                    redkey     = "HTOpca_",
                    graph      = "HTO.pca.nn",
                    red2       = "HTO.pca.umap",
                    redkey2    = "HTOpcaUMAP_",
                    dim        = 2:23,
                    res        = 0.8){
                          DefaultAssay(Seurat) = assay
                          Seurat               = ScaleData(Seurat , assay = assay)
                          Seurat               = RunPCA(Seurat       , assay=assay, reduction.name = red, 
                                                        reduction.key = redkey, verbose = FALSE, features = rownames(Seurat[[assay]]) )
                          Seurat               = FindNeighbors(Seurat, assay=assay, reduction=red,graph.name = graph, dims = dim )
                          Seurat               = FindClusters(Seurat, graph.name = graph , resolution = res )
                          Seurat               = RunUMAP(Seurat, dims = dim , reduction = red ,reduction.name = red2 , reduction.key = redkey2 )
                          Seurat
                         }

######################################################
################# Seurat HTO heatmap   ###############
######################################################

HTO_heatmap = function(Seurat        = DOGMA_pass,
                       HTO_list      = HTO_list,
                       assay         = "HTO",
                       with.Negative = TRUE,
                       low_score_zero = FALSE,
                       low_thresh    = 1.5,
                       thresh_pattern = "half",
                       TITLE         = "title",
                       OUTPUT        = "tmp.pdf"   ){
                        HTO_data = GetAssayData(Seurat, assay = assay,  slot = "data" )  %>% t() %>% as.data.frame()
                        HTO_res  = data.frame(HTO = Seurat$HTO_final)
                        
                        if(thresh_pattern == "half"){tmp_thresh = c(take_factor(HTO_list$HTO_classification,1,"_"),take_factor(HTO_list$HTO_classification,2,"_")) %>%
                                                                    table() %>% as.data.frame() %>%
                                                                    mutate(thresh = ((1 - Freq/nrow(HTO_list))*100) %>% ceiling(.)/100 ) %>%
                                                                    mutate(thresh = (1 + thresh)/2 ) %>%
                                                                    dplyr::rename(HTO = 1) %>%
                                                                    inner_join(data.frame(HTO=colnames(HTO_data)),.,by="HTO")}
                        if(thresh_pattern == "third"){tmp_thresh = c(take_factor(HTO_list$HTO_classification,1,"_"),take_factor(HTO_list$HTO_classification,2,"_")) %>%
                                                                    table() %>% as.data.frame() %>%
                                                                    mutate(thresh = ((1 - Freq/nrow(HTO_list))*100) %>% ceiling(.)/100 ) %>%
                                                                    mutate(thresh = (2 + thresh)/3 ) %>%
                                                                    dplyr::rename(HTO = 1) %>%
                                                                    inner_join(data.frame(HTO=colnames(HTO_data)),.,by="HTO")}
                        
                        HTO_data = HTO_data[,tmp_thresh$HTO]
                        
                        for(tmp_num in 1:ncol(HTO_data)){
                            HTO_data[,tmp_num] = 4*HTO_data[,tmp_num]/quantile(HTO_data[,tmp_num],tmp_thresh$thresh[tmp_num])
                        }
                        
                        if(low_score_zero){HTO_data[HTO_data<low_thresh] = 0}
                        HTO_data[HTO_data>4] = 4
                        
                        #######################################################################################
                        if(with.Negative){ HTO_list_2 = HTO_list %>% arrange(HTO_classification) %>% rbind(.,c("Negative","","gray"))}else{
                                           HTO_list_2 = HTO_list %>% arrange(HTO_classification)                                     }
                                            
                        HTO_list_2 = HTO_list_2 %>% dplyr::filter(HTO %in% unique(Seurat$HTO_final)) ######################
                        
                        tmp_list = list()
                        for(tmp_num in 1:nrow(HTO_list_2)){
                            tmp_CB_num = sum(HTO_res$HTO == HTO_list_2$HTO[tmp_num])
                            if(tmp_CB_num == 0){tmp_list[[tmp_num]] = c()}
                            if(tmp_CB_num != 0){CellBarcode_tmp     = HTO_res %>% dplyr::filter(HTO == HTO_list_2$HTO[tmp_num]) %>% rownames(.) }
                        
                            if(tmp_CB_num >1 ){tmp_list[[tmp_num]] = hclust_order(HTO_data[CellBarcode_tmp,])}
                            if(tmp_CB_num ==1){tmp_list[[tmp_num]] = CellBarcode_tmp }
                        }
                        
                        name_list = unlist(tmp_list)
                        
                        data = HTO_data[name_list,] %>% t() %>% as.data.frame()
                        
                        df = data %>% rownames_to_column("rownames") %>% gather(.,key=colnames,value=score,-rownames)
                        
                        df2 = left_join(data.frame(CellBarcode=colnames(data)) , HTO_res %>% rownames_to_column("CellBarcode"),by="CellBarcode") %>%
                               mutate(HTO = factor(HTO,levels=HTO_list_2$HTO))
                        
                        tmp_HTO_num  = rep(NA,nrow(HTO_list_2))
                        for(i in 1:nrow(HTO_list_2)){tmp_HTO_num[i]=grep(paste0("^",HTO_list_2$HTO[i],"$"),df2$HTO)%>%max();names(tmp_HTO_num)[i]=HTO_list_2$HTO[i]}
                        
                        p1  = ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                            ggrastr::geom_tile_rast()+
                            scale_fill_gradientn(colours=PurpleAndYellow())+
                            plot_theme() +
                            scale_x_discrete(limits=unique(df$colnames)) +
                            scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
                            theme(axis.text.x=element_blank(),
                                  strip.text=element_text(size=8),
                                  axis.line = element_line(colour = "black"),
                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                            labs(x="Cell Barcode",y="HTO")
                        
                        for(i in 1:nrow(HTO_list_2)){
                          p1 = p1 +
                               geom_vline(xintercept=tmp_HTO_num[i]+0.5,col="black")
                        }
                        
                        p2 = ggplot(df2, aes(x=CellBarcode,y=1,fill=HTO))+
                              ggrastr::geom_tile_rast()+
                              scale_fill_manual(values = HTO_list_2$color)+
                              scale_y_continuous(expand=c(0,0)) +
                              scale_x_discrete(limits=colnames(data)) +
                              theme(axis.title.x=element_blank(),axis.ticks=element_blank(),
                                    axis.text=element_blank(),
                                    axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
                                    plot.title = element_text(size=10,face="bold", family = "Helvetica"),
                                    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
                              labs(title=paste0(TITLE," : ",nrow(df2),"cells"),y="sgRNA")
                        
                        p_sum = p2+p1+  plot_layout(ncol = 1, heights = c(1, 9))
                        
                        pdf_2(OUTPUT ,h=5,w=15)
                         plot(p_sum)
                        dev.off()
                      }

######################################################
################# GeneActivity      ###############
######################################################

Seurat_GeneActivity = function(Seurat = data,
                               assay      = "ATAC",
                               assay_2    = "GeneActivity"){
                          DefaultAssay(Seurat) = assay
                          gene.activities      = GeneActivity(Seurat,assay=assay)
                          Seurat[[assay_2]]    = CreateAssayObject(counts = gene.activities)
                          Seurat               = NormalizeData(Seurat, assay=assay_2,normalization.method = "LogNormalize", scale.factor = 10000)
                          Seurat
                         }


######################################################
################# Seurat_summary_stat       ##########
######################################################

Seurat_pct_summary = function(Seurat     = Seurat,
                              assay      = "RNA", 
                              group.by   = "HTO_final_2",
                              sep = 1
                              ){
                                if(sep==1){df_tmp = GetAssayData(Seurat,assay=assay,slot="data") %>% 
                                                     as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("CellBarcode")}
                                if(sep>1){tmp_sep_length = ceiling(ncol(Seurat[[assay]]) / sep)
                                          df_tmp_list = list() 
                                          for(tmp_num in 1:sep){tmp_from = (tmp_num -1) * tmp_sep_length + 1
                                                                tmp_to   = min(tmp_num * tmp_sep_length,ncol(Seurat[[assay]]))
                                                                df_tmp_list[[tmp_num]] = GetAssayData(Seurat[,tmp_from:tmp_to],assay=assay,slot="data") %>% 
                                                                                          as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("CellBarcode")
                                                                }
                                          df_tmp = do.call(rbind,df_tmp_list)
                                          rm(df_tmp_list)
                                         }
                                
                                tmp_meta = Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                            dplyr::select(c("CellBarcode",group.by))
                                colnames(tmp_meta)[2] = "target"
                                
                                pct_check    = function(x){sum(x>0)/length(x)}

                                df_tmp_2 = full_join(tmp_meta , df_tmp , by = "CellBarcode") %>% 
                                            column_to_rownames("CellBarcode") %>% 
                                            group_by(target)
                                
                                pct_data = df_tmp_2 %>% summarise(across(everything(),pct_check))
                                exp_data = df_tmp_2 %>% summarise(across(everything(),mean))
                                med_data = df_tmp_2 %>% summarise(across(everything(),median))
                                len_data = df_tmp_2 %>% summarise(across(everything(),length))
                                
                                tmp_function = function(data = data , name = "function"){
                                                    data_2 = data %>%
                                                              as.data.frame(.) %>% 
                                                              pivot_longer(., col = -target, names_to = "Gene", values_to = name) %>% 
                                                              as.data.frame(.)
                                                   }
                                
                                pct_data_2 = tmp_function(data = pct_data , name = "pct") 
                                exp_data_2 = tmp_function(data = exp_data , name = "exp") 
                                med_data_2 = tmp_function(data = med_data , name = "med") 
                                len_data_2 = tmp_function(data = len_data , name = "len") 
                                
                                data_sum = do.call(function(x,y){full_join(x,y,by=c("target","Gene"))},list(pct_data_2,exp_data_2,med_data_2,len_data_2))
                                data_sum
                              }

######################################################
################# Seurat normalization ###############
######################################################
zscore_to_ctrl = function(Seurat = data,
                          assay  = "RNA",
                          slot = "data",
                          ident  = "HTO_final_2",
                          ctrl   = "NTC"
                          ){
                            Seurat_RNA   = GetAssayData(Seurat, assay = assay,   slot = slot)
                            
                            Seurat_RNA   = Seurat_RNA %>%
                                            t() %>% as.data.frame() %>% 
                                            rownames_to_column("CellBarcode") %>%
                                            inner_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                            dplyr::select("CellBarcode",ident), . , by="CellBarcode") %>%
                                            column_to_rownames("CellBarcode")
                            colnames(Seurat_RNA)[1] = "target"
                            
                            Seurat_RNA_summary  = Seurat_RNA  %>%
                                                   group_by(target) %>%
                                                   summarise(across(everything(),list(mean = mean, sd = sd) )) %>%
                                                   as.data.frame()
                            rownames(Seurat_RNA_summary) = c()
                            
                            Seurat_RNA_mean = Seurat_RNA_summary %>% 
                                               column_to_rownames("target") %>%
                                               dplyr::select(ends_with("_mean"))
                            Seurat_RNA_sd   = Seurat_RNA_summary %>% 
                                               column_to_rownames("target") %>%
                                               dplyr::select(ends_with("_sd"))
                            tmp           = Seurat_RNA_mean %>% t() %>% as.data.frame()
                            tmp_2         = tmp - tmp[,ctrl]
                            Seurat_RNA_zscore = tmp_2/Seurat_RNA_sd[ctrl,] %>% as.numeric()
                            rownames(Seurat_RNA_zscore) = gsub("_mean$","",rownames(Seurat_RNA_zscore))
                            Seurat_RNA_zscore
                           }

zscore_to_ctrl_sep = function(Seurat = data,
                          assay  = "RNA",
                          slot = "data",
                          ident  = "HTO_final_2",
                          ctrl   = "NTC",
                          tmp_from = tmp_from,
                          tmp_to = tmp_to
                          ){
                            Seurat_RNA   = GetAssayData(Seurat, assay = assay,   slot = slot)
                            
                            Seurat_RNA   = Seurat_RNA[tmp_from:tmp_to,,drop=F] %>%
                                            t() %>% as.data.frame() %>% 
                                            rownames_to_column("CellBarcode") %>%
                                            inner_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                            dplyr::select("CellBarcode",ident), . , by="CellBarcode") %>%
                                            column_to_rownames("CellBarcode")
                            colnames(Seurat_RNA)[1] = "target"
                            
                            Seurat_RNA_summary  = Seurat_RNA  %>%
                                                   group_by(target) %>%
                                                   summarise(across(everything(),list(mean = mean, sd = sd) )) %>%
                                                   as.data.frame()
                            rownames(Seurat_RNA_summary) = c()
                            
                            Seurat_RNA_mean = Seurat_RNA_summary %>% 
                                               column_to_rownames("target") %>%
                                               dplyr::select(ends_with("_mean"))
                            Seurat_RNA_sd   = Seurat_RNA_summary %>% 
                                               column_to_rownames("target") %>%
                                               dplyr::select(ends_with("_sd"))
                            tmp           = Seurat_RNA_mean %>% t() %>% as.data.frame()
                            tmp_2         = tmp - tmp[,ctrl]
                            Seurat_RNA_zscore = tmp_2/Seurat_RNA_sd[ctrl,] %>% as.numeric()
                            rownames(Seurat_RNA_zscore) = gsub("_mean$","",rownames(Seurat_RNA_zscore))
                            Seurat_RNA_zscore
                           }

zscore_slot    = function(Seurat = data,
                          assay  = "RNA",
                          slot = "data",
                          ident  = "HTO_final_2",
                          sep_num  = 1
                          ){
                            Seurat_orig   = GetAssayData(Seurat, assay = assay,   slot = slot)
                            
                            res_list = list()
                            for(tmp_num in 1:sep_num){
                                    paste0(tmp_num ,"/",sep_num," START") %>% print(.)
                                    each_num = ceiling( nrow(Seurat_orig) / sep_num )
                                    tmp_from = each_num * (tmp_num - 1) + 1
                                    tmp_to   = min(each_num * tmp_num , nrow(Seurat_orig))
    
                                    if(tmp_from<=tmp_to){
                                      Seurat_RNA = Seurat_orig[tmp_from:tmp_to,]
                                      Seurat_RNA   = Seurat_RNA %>%
                                                      t() %>% as.data.frame() %>% 
                                                      rownames_to_column("CellBarcode") %>%
                                                      inner_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                                      dplyr::select("CellBarcode",ident), . , by="CellBarcode") %>%
                                                      column_to_rownames("CellBarcode")
                                      colnames(Seurat_RNA)[1] = "target"
                                  
                                      Seurat_RNA_summary  = Seurat_RNA  %>%
                                                             group_by(target) %>%
                                                             summarise(across(everything(),list(mean = mean, sd = sd) )) %>%
                                                             as.data.frame()
                                      rownames(Seurat_RNA_summary) = c()
                                      
                                      Seurat_RNA_mean = Seurat_RNA_summary %>% 
                                                         column_to_rownames("target") %>%
                                                         dplyr::select(ends_with("_mean"))
          
                                      ########################
                                      Seurat_tmp_summary  = Seurat_RNA  %>%
                                                             dplyr::select(-target) %>%
                                                             summarise(across(everything(),list(mean = mean, sd = sd) )) %>%
                                                             as.data.frame()
                                      
                                      Seurat_tmp_mean = Seurat_tmp_summary %>% 
                                                         dplyr::select(ends_with("_mean"))
          
                                      Seurat_tmp_sd   = Seurat_tmp_summary %>% 
                                                         dplyr::select(ends_with("_sd"))
          
                                      tmp           = Seurat_tmp_mean %>% t() %>% as.data.frame()
                                      tmp_2         = Seurat_RNA_mean %>% t() %>% as.data.frame() - tmp[,1]
                                      Seurat_RNA_zscore = tmp_2/Seurat_tmp_sd[1,] %>% as.numeric()
                                      rownames(Seurat_RNA_zscore) = gsub("_mean$","",rownames(Seurat_RNA_zscore))
                                      res_list[[tmp_num]] = Seurat_RNA_zscore
                                    }
                                }
                            if(sep_num==1){res = res_list[[1]]}else{res = do.call(rbind,res_list) }
                            res
                           }

######################################################
################# Seurat summarize     ###############
######################################################
mean_assay = function(Seurat = ASAP,
                      assay  = "ADT",
                      slot   = "data",
                      ident  = "HTO_final_2",
                      with.filter = FALSE, quantile = 0.95, omit_target = NA, include_target=NA,
                      sep_num=1
                      ){
                        if(with.filter){
                          tmp_counts = GetAssayData(Seurat, assay = assay,   slot = "counts")
                          tmp_data   = GetAssayData(Seurat, assay = assay,   slot = "data"  )
                          
                          tmp_counts = tmp_counts[setdiff(rownames(tmp_counts),omit_target),]
                          tmp_data   = tmp_data[setdiff(rownames(tmp_data),omit_target),]

                          if(!is.na(include_target[1])){
                            tmp_counts = tmp_counts[intersect(rownames(tmp_counts),include_target),]
                            tmp_data   = tmp_data[intersect(rownames(tmp_data),include_target),]
                          }
                          
                          tmp_counts_2 = tmp_counts %>% t() %>% as.data.frame() %>% rownames_to_column("CellBarcode") %>% 
                                          left_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% dplyr::select("CellBarcode",ident),
                                                    . , by="CellBarcode") %>%
                                          dplyr::rename(ident = 2)
                          colnames(tmp_counts_2)[3:ncol(tmp_counts_2)] = paste0("Target_",colnames(tmp_counts_2)[3:ncol(tmp_counts_2)])
                          
                          quantile.filter = function(x){quantile(x,quantile)}
                          df_tmp   = tmp_counts_2 %>%
                                        group_by(ident) %>%
                                        summarise(across(starts_with(c("Target_")),list(quantile.filter  = quantile.filter ) )) %>%
                                        as.data.frame()
                          df_tmp_2 = df_tmp[,-1] %>% colSums()
                          LowExp_target = names(df_tmp_2)[df_tmp_2==0] %>% gsub("Target_","",.) %>% gsub("_quantile.filter","",.)
                          
                          tmp_counts = tmp_counts[setdiff(rownames(tmp_counts),LowExp_target),]
                          tmp_data   = tmp_data[setdiff(rownames(tmp_data),LowExp_target),]
                          
                          Seurat[["tmpfilt"]] = CreateAssayObject(counts = tmp_counts )
                          Seurat               = SetAssayData(object=Seurat, assay = "tmpfilt" , slot="data", 
                                                              new.data = tmp_data )
                          Seurat               = ScaleData( Seurat ,    assay = "tmpfilt")
                          assay = "tmpfilt"
                        }

                        Seurat_data  = GetAssayData(Seurat, assay = assay,   slot = slot)
                        
                        res_list = list()
                        for(tmp_num in 1:sep_num){
                            each_num = ceiling( nrow(Seurat_data) / sep_num )
                            tmp_from = each_num * (tmp_num - 1) + 1
                            tmp_to   = min(each_num * tmp_num , nrow(Seurat_data))

                            if(tmp_from<=tmp_to){
                               Seurat_data_p  = Seurat_data[tmp_from:tmp_to,] %>%
                                                  t() %>% as.data.frame() %>% 
                                                  rownames_to_column("CellBarcode") %>%
                                                  inner_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                                  dplyr::select("CellBarcode",ident), . , by="CellBarcode") %>%
                                                  column_to_rownames("CellBarcode")
                               colnames(Seurat_data_p)[1] = "target"
                            
                               Seurat_data_summary  = Seurat_data_p  %>%
                                                      group_by(target) %>%
                                                      summarise(across(everything(),list(mean = mean) )) %>%
                                                      as.data.frame()
                               rownames(Seurat_data_summary) = c()
                            
                               Seurat_data_mean = Seurat_data_summary %>% 
                                                  column_to_rownames("target") %>%
                                                  dplyr::select(ends_with("_mean"))
                               colnames(Seurat_data_mean) = gsub("_mean$","",colnames(Seurat_data_mean))
                               res_list[[tmp_num]] = Seurat_data_mean
                            }
                            show.progress(tmp_num, 1:sep_num)
                        }
                        if(sep_num==1){res = res_list[[1]]}else{res = do.call(cbind,res_list) }
                        res
                      }

######################################################
################# Seurat summarize     ###############
######################################################
median_assay = function(Seurat = ASAP,
                        assay  = "ADT",
                        slot   = "data",
                        ident  = "HTO_final_2",
                        with.filter = FALSE, quantile = 0.95, omit_target = NA, include_target=NA,
                        sep_num=1
                        ){
                          if(with.filter){
                            tmp_counts = GetAssayData(Seurat, assay = assay,   slot = "counts")
                            tmp_data   = GetAssayData(Seurat, assay = assay,   slot = "data"  )
                            
                            tmp_counts = tmp_counts[setdiff(rownames(tmp_counts),omit_target),]
                            tmp_data   = tmp_data[setdiff(rownames(tmp_data),omit_target),]
  
                            if(!is.na(include_target[1])){
                              tmp_counts = tmp_counts[intersect(rownames(tmp_counts),include_target),]
                              tmp_data   = tmp_data[intersect(rownames(tmp_data),include_target),]
                            }
                            
                            tmp_counts_2 = tmp_counts %>% t() %>% as.data.frame() %>% rownames_to_column("CellBarcode") %>% 
                                            left_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% dplyr::select("CellBarcode",ident),
                                                      . , by="CellBarcode") %>%
                                            dplyr::rename(ident = 2)
                            colnames(tmp_counts_2)[3:ncol(tmp_counts_2)] = paste0("Target_",colnames(tmp_counts_2)[3:ncol(tmp_counts_2)])
                            
                            quantile.filter = function(x){quantile(x,quantile)}
                            df_tmp   = tmp_counts_2 %>%
                                          group_by(ident) %>%
                                          summarise(across(starts_with(c("Target_")),list(quantile.filter  = quantile.filter ) )) %>%
                                          as.data.frame()
                            df_tmp_2 = df_tmp[,-1] %>% colSums()
                            LowExp_target = names(df_tmp_2)[df_tmp_2==0] %>% gsub("Target_","",.) %>% gsub("_quantile.filter","",.)
                            
                            tmp_counts = tmp_counts[setdiff(rownames(tmp_counts),LowExp_target),]
                            tmp_data   = tmp_data[setdiff(rownames(tmp_data),LowExp_target),]
                            
                            Seurat[["tmpfilt"]] = CreateAssayObject(counts = tmp_counts )
                            Seurat               = SetAssayData(object=Seurat, assay = "tmpfilt" , slot="data", 
                                                                new.data = tmp_data )
                            Seurat               = ScaleData( Seurat ,    assay = "tmpfilt")
                            assay = "tmpfilt"
                          }

                          Seurat_data  = GetAssayData(Seurat, assay = assay,   slot = slot)
                          
                          res_list = list()
                          for(tmp_num in 1:sep_num){
                              each_num = ceiling( nrow(Seurat_data) / sep_num )
                              tmp_from = each_num * (tmp_num - 1) + 1
                              tmp_to   = min(each_num * tmp_num , nrow(Seurat_data))
  
                              if(tmp_from<=tmp_to){
                                Seurat_data_p  = Seurat_data[tmp_from:tmp_to,] %>%
                                                t() %>% as.data.frame() %>% 
                                                rownames_to_column("CellBarcode") %>%
                                                inner_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                                dplyr::select("CellBarcode",ident), . , by="CellBarcode") %>%
                                                column_to_rownames("CellBarcode")
                                colnames(Seurat_data_p)[1] = "target"
                                
                                Seurat_data_summary  = Seurat_data_p  %>%
                                                       group_by(target) %>%
                                                       summarise(across(everything(),list(median = median) )) %>%
                                                       as.data.frame()
                                rownames(Seurat_data_summary) = c()
                                
                                Seurat_data_mean = Seurat_data_summary %>% 
                                                   column_to_rownames("target") %>%
                                                   dplyr::select(ends_with("_median"))
                                colnames(Seurat_data_mean) = gsub("_median$","",colnames(Seurat_data_mean))
                                res_list[[tmp_num]] = Seurat_data_mean
                             }
                            show.progress(tmp_num, 1:sep_num)
                          }
                          if(sep_num==1){res = res_list[[1]]}else{res = do.call(cbind,res_list) }
                          res
                         }

######################################################
################# Seurat summarize     ###############
######################################################
zscore_assay = function(Seurat = ASAP,
                        assay  = "ADT",
                        slot   = "data",
                        ident  = "HTO_final_2",
                        with.filter = FALSE, quantile = 0.95, omit_target = NA, include_target=NA,
                        sep_num=1
                        ){
                        if(with.filter){
                          tmp_counts = GetAssayData(Seurat, assay = assay,   slot = "counts")
                          tmp_data   = GetAssayData(Seurat, assay = assay,   slot = "data"  )
                          
                          tmp_counts = tmp_counts[setdiff(rownames(tmp_counts),omit_target),]
                          tmp_data   = tmp_data[setdiff(rownames(tmp_data),omit_target),]

                          if(!is.na(include_target[1])){
                            tmp_counts = tmp_counts[intersect(rownames(tmp_counts),include_target),]
                            tmp_data   = tmp_data[intersect(rownames(tmp_data),include_target),]
                          }
                          
                          tmp_counts_2 = tmp_counts %>% t() %>% as.data.frame() %>% rownames_to_column("CellBarcode") %>% 
                                          left_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% dplyr::select("CellBarcode",ident),
                                                    . , by="CellBarcode") %>%
                                          dplyr::rename(ident = 2)
                          colnames(tmp_counts_2)[3:ncol(tmp_counts_2)] = paste0("Target_",colnames(tmp_counts_2)[3:ncol(tmp_counts_2)])
                          
                          quantile.filter = function(x){quantile(x,quantile)}
                          df_tmp   = tmp_counts_2 %>%
                                        group_by(ident) %>%
                                        summarise(across(starts_with(c("Target_")),list(quantile.filter  = quantile.filter ) )) %>%
                                        as.data.frame()
                          df_tmp_2 = df_tmp[,-1] %>% colSums()
                          LowExp_target = names(df_tmp_2)[df_tmp_2==0] %>% gsub("Target_","",.) %>% gsub("_quantile.filter","",.)
                          
                          tmp_counts = tmp_counts[setdiff(rownames(tmp_counts),LowExp_target),]
                          tmp_data   = tmp_data[setdiff(rownames(tmp_data),LowExp_target),]
                          
                          Seurat[["tmpfilt"]] = CreateAssayObject(counts = tmp_counts )
                          Seurat               = SetAssayData(object=Seurat, assay = "tmpfilt" , slot="data", 
                                                              new.data = tmp_data )
                          Seurat               = ScaleData( Seurat ,    assay = "tmpfilt")
                          assay = "tmpfilt"
                        }

                        Seurat_data  = GetAssayData(Seurat, assay = assay,   slot = slot)
                        
                        res_list = list()
                        for(tmp_num in 1:sep_num){
                            each_num = ceiling( nrow(Seurat_data) / sep_num )
                            tmp_from = each_num * (tmp_num - 1) + 1
                            tmp_to   = min(each_num * tmp_num , nrow(Seurat_data))

                            if(tmp_from<=tmp_to){
                               Seurat_data_p  = Seurat_data[tmp_from:tmp_to,] %>%
                                                  t() %>% as.data.frame() %>% 
                                                  rownames_to_column("CellBarcode") %>%
                                                  inner_join(Seurat@meta.data %>% rownames_to_column("CellBarcode") %>% 
                                                  dplyr::select("CellBarcode",ident), . , by="CellBarcode") %>%
                                                  column_to_rownames("CellBarcode")
                               colnames(Seurat_data_p)[1] = "target"
                            
                               Seurat_data_summary  = Seurat_data_p  %>%
                                                      group_by(target) %>%
                                                      summarise(across(everything(),list(mean = mean) )) %>%
                                                      as.data.frame()
                               rownames(Seurat_data_summary) = c()

                               Seurat_data_mean = Seurat_data_summary %>% 
                                                  column_to_rownames("target") %>%
                                                  dplyr::select(ends_with("_mean"))
                               colnames(Seurat_data_mean) = gsub("_mean$","",colnames(Seurat_data_mean))
                               ##########################
                               Seurat_tmp_summary  = Seurat_data_p %>%
                                                      dplyr::select(-target) %>%
                                                      summarise(across(everything(),list(mean = mean, sd = sd) )) %>%
                                                      as.data.frame()
                               Seurat_tmp_mean = Seurat_tmp_summary %>% 
                                                  dplyr::select(ends_with("_mean")) %>%
                                                  t(.) %>%
                                                  as.data.frame(.)
                               Seurat_tmp_sd   = Seurat_tmp_summary %>% 
                                                  dplyr::select(ends_with("_sd")) %>%
                                                  t(.) %>%
                                                  as.data.frame(.)
                               rownames(Seurat_tmp_mean) = rownames(Seurat_tmp_mean) %>% gsub("_mean","",.)
                               rownames(Seurat_tmp_sd)   = rownames(Seurat_tmp_sd)   %>% gsub("_sd",  "",.)
                               colnames(Seurat_tmp_mean) = "mean"
                               colnames(Seurat_tmp_sd)   = "sd"
                               Seurat_tmp = cbind(Seurat_tmp_mean,Seurat_tmp_sd)

                               Seurat_tmp_2      = Seurat_data_mean %>% t(.) %>% as.data.frame(.) - Seurat_tmp[,"mean"]
                               Seurat_data_zscore = Seurat_tmp_2/Seurat_tmp[,"sd"]
                               Seurat_data_zscore = Seurat_data_zscore %>% t(.) %>% as.data.frame(.)
                               res_list[[tmp_num]] = Seurat_data_zscore
                            }
                            show.progress(tmp_num, 1:sep_num)
                        }
                        if(sep_num==1){res = res_list[[1]]}else{res = do.call(cbind,res_list) }
                        res
                      }

######################################################
################# DEG by edgeR           #############
######################################################

edgeRQLF_Seurat  =function(Seurat = CITE,
                           group  = "HTO_final_2",
                           ctrl   = "NTC",
                           assay  = "RNA",
                           model_batch  = NA,
                           omit_gene = NA,
                           OUTPUT = "output.txt"
                           ){
                            library("SingleCellExperiment") %>% suppressMessages()
                            library("edgeR")                %>% suppressMessages()
                            library("scran")                %>% suppressMessages()
                            
                            # count matrix 
                            counts <- GetAssayData(Seurat , slot="counts" , assay = assay) %>% as.matrix(.)
                            counts <- counts[rowSums(counts) >= 1, ]
                            if(!is.na(omit_gene[1])){counts = counts[setdiff(rownames(counts),omit_gene)%>%sort(.),]}
                            
                            # subset the meta data fro filtered gene/cells
                            if(is.na(model_batch)){metadata           = Seurat@meta.data[,group,drop=F]
                                                   colnames(metadata) = c("group")
                                             }else{metadata           = Seurat@meta.data[,c(group,model_batch),drop=F]
                                                   colnames(metadata) = c("group","model_batch")
                                           }
                            metadata <- metadata[colnames(counts),,drop=F]
                            
                            ### Make single cell experiment
                            sce <- SingleCellExperiment(assays = counts, 
                                                        colData = metadata)
                            
                            ## convert to edgR object
                            dge      = convertTo(sce, type="edgeR", assay.type = 1)
                            meta_dge = dge$samples %>%
                                        dplyr::select(c("lib.size","norm.factors")) %>%
                                        cbind(.,metadata) %>%
                                        mutate(group = factor(group) %>% relevel(.,ref=ctrl))
                            dge$samples = meta_dge
                            
                            dge    = calcNormFactors(dge)
                            if(is.na(model_batch)){design = model.matrix(~0+group,             data=dge$samples)
                                             }else{design = model.matrix(~0+group+model_batch, data=dge$samples) }
                            y      = estimateGLMCommonDisp ( dge, design )
                            y      = estimateGLMTrendedDisp( y  , design )
                            y      = estimateGLMTagwiseDisp( y  , design )
                            fit    = glmQLFit(y, design = design)
                            if(is.na(model_batch)){qlf    = glmQLFTest(fit, contrast = c(-1, 1))
                                             }else{qlf    = glmQLFTest(fit, contrast = c(-1, 1, 0)) }
                            tt     = topTags(qlf, n = Inf) %>% as.data.frame(.)
                            write.table_n_2(tt,"Gene",OUTPUT)
                           }

edgeRQLF_f_Seurat=function(Seurat = CITE,
                           group  = "HTO_final_2",
                           ctrl   = "NTC",
                           OUTPUT = "output.txt"
                           ){
                            library("SingleCellExperiment") %>% suppressMessages()
                            library("edgeR")                %>% suppressMessages()
                            library("scran")                %>% suppressMessages()
                            
                            # count matrix 
                            counts <- as.matrix(Seurat@assays$RNA@counts)
                            counts <- counts[rowSums(counts) >= 1, ]
                            
                            # subset the meta data fro filtered gene/cells
                            metadata <- Seurat@meta.data[,group,drop=F]
                            metadata <- metadata[colnames(counts),,drop=F]
                            colnames(metadata) = c("group")
                            
                            ### Make single cell experiment
                            sce <- SingleCellExperiment(assays = counts, 
                                                        colData = metadata)
                            
                            ## convert to edgR object
                            dge      = convertTo(sce, type="edgeR", assay.type = 1)
                            meta_dge = dge$samples %>%
                                        dplyr::select(c("lib.size","norm.factors")) %>%
                                        cbind(.,metadata) %>%
                                        mutate(group = factor(group) %>% relevel(.,ref=ctrl))
                            dge$samples = meta_dge
                            
                            dge    = calcNormFactors(dge)
                            design = model.matrix(~0+group, data=dge$samples)
                            y      = estimateGLMCommonDisp ( dge, design )
                            y      = estimateGLMTrendedDisp( y  , design )
                            y      = estimateGLMTagwiseDisp( y  , design )
                            fit    = glmQLFit(y, design = design)
                            qlf    = glmQLFTest(fit, contrast = c(-1, 1))
                            tt     = topTags(qlf, n = Inf) %>% as.data.frame(.)
                            write.table_n_2(tt,"Gene",OUTPUT)
                           }

edgeRQLF_Seurat_n  =function(Seurat = CITE,
                           group  = "HTO_final_2",
                           ctrl   = "CTRL",
                           trgt   = "TRGT",
                           assay  = "RNA",
                           model_batch  = NA,
                           omit_gene = NA,
                           OUTPUT = "output.txt"
                           ){
                            library("SingleCellExperiment") %>% suppressMessages()
                            library("edgeR")                %>% suppressMessages()
                            library("scran")                %>% suppressMessages()
                            
                            # count matrix 
                            counts <- GetAssayData(Seurat , layer="counts" , assay = assay) %>% as.matrix(.)
                            counts <- counts[rowSums(counts) >= 1, ]
                            if(!is.na(omit_gene[1])){counts = counts[setdiff(rownames(counts),omit_gene)%>%sort(.),]}
                            
                            # subset the meta data fro filtered gene/cells
                            if(is.na(model_batch)){metadata           = Seurat@meta.data[,group,drop=F]
                                                   colnames(metadata) = c("group")
                                             }else{metadata           = Seurat@meta.data[,c(group,model_batch),drop=F]
                                                   colnames(metadata) = c("group","model_batch")
                                           }
                            metadata <- metadata[colnames(counts),,drop=F]
                            
                            ### Make single cell experiment
                            sce <- SingleCellExperiment(assays = counts, 
                                                        colData = metadata)
                            
                            ## convert to edgR object
                            dge      = convertTo(sce, type="edgeR", assay.type = 1)
                            meta_dge = dge$samples %>%
                                        dplyr::select(c("lib.size","norm.factors")) %>%
                                        cbind(.,metadata) ## %>%
                                        ## mutate(group = factor(group) %>% relevel(.,ref=ctrl))
                            dge$samples = meta_dge
                            
                            dge    = calcNormFactors(dge)
                            if(is.na(model_batch)){design = model.matrix(~0+group,             data=dge$samples)
                                             }else{design = model.matrix(~0+group+model_batch, data=dge$samples) }

                            contrast_orig = rep(0, ncol(design)) %>%
                                            {.[colnames(design) %in% paste0("group",ctrl)]=-1;
                                             .[colnames(design) %in% paste0("group",trgt)]=1;. } %>%
                                            {names(.)=colnames(design);.}
                            tmp_sumstat  = table(metadata$group) %>% as.data.frame() %>%
                                            mutate(Var1 = paste0("group",Var1)) %>%
                                            left_join(data.frame(Var1=names(contrast_orig), ctrl_trgt = contrast_orig), ., by="Var1") %>% 
                                            mutate(Freq = case_when( is.na(Freq)~0, TRUE~Freq) )  %>% 
                                            mutate(Freq_2 = Freq * ctrl_trgt)
                            nega_sum = sum(tmp_sumstat$Freq_2[tmp_sumstat$Freq_2<0])
                            posi_sum = sum(tmp_sumstat$Freq_2[tmp_sumstat$Freq_2>0])
                            contrast_tmp  = tmp_sumstat %>%
                                             mutate(contrast =case_when( (ctrl_trgt>0)~Freq/posi_sum, (ctrl_trgt<0)~Freq/nega_sum, TRUE~0 ) ) %>%
                                             pull(contrast)

                            y      = estimateGLMCommonDisp ( dge, design )
                            y      = estimateGLMTrendedDisp( y  , design )
                            y      = estimateGLMTagwiseDisp( y  , design )
                            fit    = glmQLFit(y, design = design)
                            qlf    = glmQLFTest(fit, contrast = contrast_tmp)
                            tt     = topTags(qlf, n = Inf) %>% as.data.frame(.)
                            write.table_n_2(tt,"Gene",OUTPUT)
                           }

edgeRQLF_f_ASAP =function(Seurat = ASAP,
                           group  = "HTO_final_2",
                           ctrl   = "NTC",
                           OUTPUT = "output.txt"
                           ){
                            library("SingleCellExperiment") %>% suppressMessages()
                            library("edgeR")                %>% suppressMessages()
                            library("scran")                %>% suppressMessages()
                            
                            # count matrix 
                            counts <- as.matrix(Seurat@assays$GeneScore@counts)
                            counts <- counts[rowSums(counts) > 0, ]
                            
                            # subset the meta data fro filtered gene/cells
                            metadata <- Seurat@meta.data[,group,drop=F]
                            metadata <- metadata[colnames(counts),,drop=F]
                            colnames(metadata) = c("group")
                            
                            ### Make single cell experiment
                            sce <- SingleCellExperiment(assays = counts, 
                                                        colData = metadata)
                            
                            ## convert to edgR object
                            dge      = convertTo(sce, type="edgeR", assay.type = 1)
                            meta_dge = dge$samples %>%
                                        dplyr::select(c("lib.size","norm.factors")) %>%
                                        cbind(.,metadata) %>%
                                        mutate(group = factor(group) %>% relevel(.,ref=ctrl))
                            dge$samples = meta_dge
                            
                            dge    = calcNormFactors(dge)
                            design = model.matrix(~0+group, data=dge$samples)
                            y      = estimateGLMCommonDisp ( dge, design )
                            y      = estimateGLMTrendedDisp( y  , design )
                            y      = estimateGLMTagwiseDisp( y  , design )
                            fit    = glmQLFit(y, design = design)
                            qlf    = glmQLFTest(fit, contrast = c(-1, 1))
                            tt     = topTags(qlf, n = Inf) %>% as.data.frame(.)
                            write.table_n_2(tt,"Gene",OUTPUT)
                           }

######################################################
################# PathwayEnrichment_hg38 #############
######################################################

PathwayEnrichment_hg38 = function(gene_list = gene_name,
                                  full_gene_list = NA,
                                  name,
                                  OUTPUT = "Pathway_Enrichment",
                                  simplify = FALSE){
   suppressMessages(library(clusterProfiler))
   dir.create_p(OUTPUT)

   gene_ENTREZ_list = fread_FT("~/reference_home/enrichment_analysis/hg38_refGene_gene_w_ENTREZ_ID_list.txt")
   gene_ENTREZ_list = gene_ENTREZ_list[!is.na(gene_ENTREZ_list$ENTREZID),]
   if(!is.na(full_gene_list[1])){gene_ENTREZ_list = gene_ENTREZ_list %>% dplyr::filter(SYMBOL %in% full_gene_list)}
   gene_list        = gene_list[is.element(gene_list,gene_ENTREZ_list$SYMBOL)]
   
   #################################################

   if(length(gene_list)>0){
       target_tmp = inner_join(gene_ENTREZ_list,data.frame(SYMBOL=gene_list),by="SYMBOL")
       enrichGO  = enrichGO(gene     = target_tmp$ENTREZID         %>% as.character(),
                            universe = gene_ENTREZ_list$ENTREZID   %>% as.character(),
                            OrgDb    = "org.Hs.eg.db",
                            ont      = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                            minGSSize = 10, maxGSSize = 500,
                            readable = TRUE)
       write.table_FT_2(enrichGO,paste0(OUTPUT,"/",name,"_enrichGO.txt"))
       enrichGO.simple = simplify(enrichGO)
       write.table_FT_2(enrichGO.simple,paste0(OUTPUT,"/",name,"_enrichGO_simple.txt"))
      
      if(simplify){enrichGO = enrichGO.simple}
      if(nrow(enrichGO)>0){
          p_bar = barplot(enrichGO, showCategory=20)+ plot_annotation(title = name)
          paste0(OUTPUT,"/",name,"_enrichGO_barplot.pdf")     %>% pdf_2(.,h=10,w=10)
           plot(p_bar)
          dev.off()
          
          p_dot = dotplot(enrichGO, showCategory=20)+ plot_annotation(title = name)
          paste0(OUTPUT,"/",name,"_enrichGO_dotplot.pdf")     %>% pdf_2(.,h=10,w=10)
           plot(p_dot)
          dev.off()
          
          paste0(OUTPUT,"/",name,"_enrichGO_plotGOgraph.pdf") %>% pdf_2(.,h=10,w=10)
           plotGOgraph(enrichGO)
          dev.off()

          ## clusterProfiler::emapplot(enrichGO)
          ## clusterProfiler::cnetplot(enrichGO, categorySize="pvalue", foldChange=geneList)
          ## goplot(ego_result.simple)
      }
   }
}

######################################################
################# PathwayEnrichment_mm10 #############
######################################################

PathwayEnrichment_mm10 = function(gene_list = gene_name,
                                  full_gene_list = NA,
                                  name,
                                  OUTPUT = "Pathway_Enrichment",
                                  simplify = FALSE){
   suppressMessages(library(clusterProfiler))
   dir.create_p(OUTPUT)

   gene_ENTREZ_list = fread_FT("~/reference_home/enrichment_analysis/mm10_refGene_gene_w_ENTREZ_ID_list.txt")
   gene_ENTREZ_list = gene_ENTREZ_list[!is.na(gene_ENTREZ_list$ENTREZID),]
   if(!is.na(full_gene_list[1])){gene_ENTREZ_list = gene_ENTREZ_list %>% dplyr::filter(SYMBOL %in% full_gene_list)}
   gene_list        = gene_list[is.element(gene_list,gene_ENTREZ_list$SYMBOL)]
   
   #################################################

   if(length(gene_list)>0){
       target_tmp = inner_join(gene_ENTREZ_list,data.frame(SYMBOL=gene_list),by="SYMBOL")
       enrichGO  = enrichGO(gene     = target_tmp$ENTREZID %>% as.character(),
                            universe = gene_ENTREZ_list$ENTREZID %>% as.character(),
                            OrgDb    = "org.Mm.eg.db",
                            ont      = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                            minGSSize = 10, maxGSSize = 500,
                            readable = TRUE)
       write.table_FT_2(enrichGO,paste0(OUTPUT,"/",name,"_enrichGO.txt"))
       enrichGO.simple = simplify(enrichGO)
       write.table_FT_2(enrichGO.simple,paste0(OUTPUT,"/",name,"_enrichGO_simple.txt"))
      
      if(simplify){enrichGO = enrichGO.simple}
      if(nrow(enrichGO)>0){
          p_bar = barplot(enrichGO, showCategory=20)+ plot_annotation(title = name)
          paste0(OUTPUT,"/",name,"_enrichGO_barplot.pdf")     %>% pdf_2(.,h=10,w=10)
           plot(p_bar)
          dev.off()
          
          p_dot = dotplot(enrichGO, showCategory=20)+ plot_annotation(title = name)
          paste0(OUTPUT,"/",name,"_enrichGO_dotplot.pdf")     %>% pdf_2(.,h=10,w=10)
           plot(p_dot)
          dev.off()
          
          paste0(OUTPUT,"/",name,"_enrichGO_plotGOgraph.pdf") %>% pdf_2(.,h=10,w=10)
           plotGOgraph(enrichGO)
          dev.off()

          ## clusterProfiler::emapplot(enrichGO)
          ## clusterProfiler::cnetplot(enrichGO, categorySize="pvalue", foldChange=geneList)
          ## goplot(ego_result.simple)
      }
   }
}

######################################################
################# PathwayEnrichment      #############
######################################################

Pathway_Enrichment_dhyper = function(OUTPUT = "cluster_enrichment",name = "NA" , target_list, full_gene_list = NA, pathway = "BioPlanet"){
 if(pathway == "IPA"){path_of_pathway = "~/reference_home/pathway/180530_pathway_gene_list_2.txt"}
 if(pathway == "BioPlanet"){path_of_pathway = "~/reference_home/BioPlanet/BioPlanet_geneset.txt"}
 PATHWAY_LIST=read.table(path_of_pathway,header=T,row.names=NULL,stringsAsFactors=F,sep="\t")

 if(is.na(full_gene_list[1])){full_gene_list = fread_FT("~/data_home/TCR100_iTreg100/data/220621_target_gene_full_list.txt") %>% pull(Gene)}

 FILTER_PATH=function(x){as.numeric_f(x[2])*as.numeric_f(x[3])!=0}
 DHYPER=function(x){1-phyper(as.numeric_f(x[2])-1, as.numeric_f(x[3]), as.numeric_f(x[4]), as.numeric_f(x[6]))}
 ODDS=function(x){as.numeric_f(x[2])*(as.numeric_f(x[3])+as.numeric_f(x[4])-as.numeric_f(x[6]))/(as.numeric_f(x[6])*(as.numeric_f(x[3])-as.numeric_f(x[2])))}

 pathway_gene_list = strsplit(PATHWAY_LIST$gene,",")

 RES=data.frame(PATHWAY=PATHWAY_LIST$pathway) %>%
      mutate(intersect_target_PATH  = NA,
             intersect_GENE_PATH =  NA,
             setdiff_GENE_PATH = NA,
             intersect_target_name = NA)

 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_target_PATH[kkk]=length(intersect(target_list,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_GENE_PATH[kkk]=length(intersect(full_gene_list,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$setdiff_GENE_PATH[kkk]=length(setdiff(full_gene_list,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_target_name[kkk]=intersect(target_list,pathway_gene_list[[kkk]]) %>% paste0(.,collapse=",")}
 RES$num_target=length(target_list)

 RES_2=RES[apply(RES,1,FILTER_PATH),]
 RES_2$Pvalue_hyper=apply(RES_2,1,DHYPER)
 RES_2$logPvalue=round(-log10(RES_2$Pvalue_hyper),digit=3)
 RES_2$Qvalue_BH=p.adjust(RES_2$Pvalue_hyper, method = "BH")
 RES_2$ODDS=apply(RES_2,1,ODDS)
 RES_2$RATIO=round(RES_2$intersect_target_PATH/RES_2$intersect_GENE_PATH,digits=3)
 RES_3=RES_2[order(RES_2$Pvalue_hyper),]
 RES_3 = RES_3[,c("PATHWAY","intersect_target_PATH","intersect_GENE_PATH",
                  "setdiff_GENE_PATH","num_target",
                  "Pvalue_hyper","logPvalue","Qvalue_BH",
                  "ODDS","RATIO","intersect_target_name")]

 if(pathway == "BioPlanet"){
    corresp = fread_FT("~/reference_home/BioPlanet/BioPlanet_PATHWAY_ID_NAME.txt")
    corresp_list = corresp$PATHWAY_NAME
    names(corresp_list) = corresp$PATHWAY_ID
    RES_3 = RES_3 %>% 
             mutate(PATHWAY = PATHWAY %>% str_replace_all(corresp_list) )      }

 write.table_FT_2(RES_3,paste0(OUTPUT,"/full_res/",name,"_Pathway_enrichment_full.txt"))
 write.table_FT_2(RES_3[(RES_3$Qvalue_BH<0.05)&(RES_3$intersect_GENE_PATH>5),c(1,6:11)],paste0(OUTPUT,"/res_part/",name,"_Pathway_enrichment_filter.txt"))
 }

######################################################
################# Get module score       #############
######################################################

SignatureScore = function(Seurat = Seurat, 
                           assay  = "RNA.SCT.CC",
                           OUTPUT = "SignatureScore_summary.txt",
                           gene_list_path = "/home/yutake/reference_home/enrichment_analysis/MSigDB_gene_list.rds",
                           gene_list = NULL,
                           tmp_from = NA , tmp_to = NA, tmp_seed = 1234, nbin = 24){
                     ### #****************************************#
                     ### # Prep of signature score
                     ### #****************************************#
                     ### 
                     ### HALLMARK = read.gmt("~/reference/MSigDB/h.all.v7.4.symbols.gmt")
                     ### REACTOME = read.gmt("~/reference/MSigDB/c2.cp.reactome.v7.4.symbols.gmt")
                     ### KEGG     = read.gmt("~/reference/MSigDB/c2.cp.kegg.v7.4.symbols.gmt")
                     ### WikiPath = read.gmt("~/reference/MSigDB/c2.cp.wikipathways.v7.4.symbols.gmt")
                     ### GO_BP    = read.gmt("~/reference/MSigDB/c5.go.bp.v7.4.symbols.gmt")
                     ### 
                     ### geneset      = rbind(HALLMARK,REACTOME,KEGG,WikiPath,GO_BP)
                     ### geneset$term = as.character(geneset$term)
                     ### gene_set_name_list = unique(geneset$term)
                     ### 
                     ### gene_list = list()
                     ### for(tmp_num in 1:length(gene_set_name_list) ){
                     ###   gene_set_tmp = gene_set_name_list[tmp_num]
                     ###   gene_list[[tmp_num]] = geneset %>% filter(term==gene_set_tmp) %>%
                     ###                           pull(gene) %>% 
                     ###                           as.character() %>% sort()
                     ### }
                     ### 
                     ### names(gene_list) = gene_set_name_list
                     ### saveRDS(gene_list, file= paste0("/home/yutake/reference_home/enrichment_analysis/MSigDB_gene_list.rds"))
                     
                     ################
                     if(!is.na(gene_list_path)){gene_list = readRDS(gene_list_path)}
                     
                     if(is.na(tmp_from)){tmp_from  = 1}
                     if(is.na(tmp_to)){tmp_to      = length(gene_list) }
                     
                     gene_list_2        = sapply(gene_list[tmp_from:tmp_to],function(x){intersect(x,rownames(Seurat[[assay]]))})
                     gene_list_3        = gene_list_2[sapply(gene_list_2,length)>=5]
                     gene_set_name_list = names(gene_list_3)
                     
                     Seurat = AddModuleScore(object = Seurat,
                                             assay = assay, 
                                             features = gene_list_3,
                                             nbin = nbin,
                                             ctrl = 100,
                                             name = "SignatureScore",
                                             seed = tmp_seed   )
                     tmp_data = Seurat@meta.data %>% dplyr::select(starts_with("SignatureScore"))
                     colnames(tmp_data) =  gene_set_name_list
                     
                     write.table_n_2(tmp_data,"CellBarcode",OUTPUT)
}

SignatureScore_single = function(Seurat = Seurat, 
                           assay  = "RNA.SCT.CC",
                           OUTPUT = "SignatureScore_summary.txt",
                           gene_list_path = "/home/yutake/reference_home/enrichment_analysis/MSigDB_gene_list.rds",
                           gene_list = NULL,
                           tmp_from = NA, tmp_seed = 1234){                     
                     ################
                     if(!is.na(gene_list_path)){gene_list = readRDS(gene_list_path)}
                     
                     if(is.na(tmp_from)){tmp_from  = 1}
                     
                     gene_list_2        = sapply(gene_list[tmp_from],function(x){intersect(x,rownames(Seurat[[assay]]))})
                     ## NULLだとerror!!!!!!!!!!!!!!!!!
                     if(length(gene_list_2)>=5){
                        gene_list_3        = gene_list_2[,1]
                        gene_set_name_list = colnames(gene_list_2)
                        
                        Seurat = AddModuleScore(object = Seurat,
                                                assay  = assay, 
                                                features = list(gene_list_3),
                                                nbin = 24,
                                                ctrl = 100,
                                                name = "SignatureScore",
                                                seed = tmp_seed   )
                        tmp_data = Seurat@meta.data %>% dplyr::select(starts_with("SignatureScore"))
                        colnames(tmp_data) =  gene_set_name_list
                        
                        write.table_n_2(tmp_data,"CellBarcode",OUTPUT)
                     }
}

######################################################
################# Calculation complexity #############
######################################################

# Function translated from java version: https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java
# Not vectorized!!! 
estimateLibrarySize <- function(nTotal, nUnique){f <- function(x, c, n) { return(c / x - 1 + exp(-n / x)) }
                                                 
                                                 m = 1
                                                 M = 100
                                                 
                                                 nDuplicates <- (nTotal - nUnique) + 1 # charity to handle only unique reads observed
                                                 
                                                 if (nUnique > nTotal | (f(m * nUnique, nUnique, nTotal) < 0) | nUnique < 0 | nTotal < 0 | nDuplicates < 0) {
                                                   message("Library size returns 0 -- invalid inputs; check this cell more closely")
                                                   return(0)
                                                 }
                                                 
                                                 while (f(M * nUnique, nUnique, nTotal) > 0) {
                                                   M <- M*10.0
                                                 }
                                                 
                                                 for(i in seq(0, 40)) {
                                                   r <-  (m + M) / 2.0
                                                   u <- f(r * nUnique, nUnique, nTotal);
                                                   if (u == 0) {
                                                     break
                                                   } else if (u > 0) {
                                                     m = r
                                                   } else if (u < 0) {
                                                     M = r
                                                   }
                                                 }
                                                 
                                                 return(round(nUnique * (m + M) / 2.0))
                                                }

######################################################
################# Jaccard Index          #############
######################################################

Jaccard = function (x, y) { return ( (intersect(x,y) %>% length()) / (union(x,y) %>% length()) ) }

JaccardIndex = function(df = df){m = matrix(data = NA, nrow = length(df), ncol = length(df))
                                 for (r in 1:length(df)) {
                                     for (c in 1:length(df)) {
                                         if (c == r) {
                                             m[r,c] = 1
                                             paste0(r,":",c) %>% print()
                                         } else if (c > r) {
                                             m[r,c] = Jaccard(df[[r]], df[[c]] )
                                             paste0(r,":",c) %>% print()
                                         } else if (r > c) {
                                             m[r,c] = m[c,r]
                                             paste0(r,":",c) %>% print()
                                         }
                                     }
                                 }
                                 m = as.data.frame(m)
                                 colnames(m) = rownames(m) = colnames(df)
                                 order = hclust_col_order(m,METHOD="complete")
                                 m[order,order] %>% as.data.frame()
                                }

######################################################
################# plot pheatmap          #############
######################################################

plot_pheatmap = function(data = data,
                         row_column = FALSE, col_column = FALSE,
                         row_annot = NA, col_annot = NA,
                         annot_palette = NA,
                         OUTPUT = "tmp.pdf",
                         TITLE  = NA,
                         heatmap_col = heatmap_col,
                         cluster_cols = TRUE, cluster_rows = TRUE,
                         cellwidth  = 10, cellheight  = 10, 
                         border_color = "black", 
                         display_numbers = FALSE, fontsize_number = 5,
                         fontsize = 10 , fontsize_row = NA, fontsize_col = NA,
                         gaps_row_num = NULL, gaps_col_num = NULL,
                         cutree_rows = NA, cutree_cols = NA,
                         clustering_method = "ward.D2",
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean"
                         ){
                          if(row_column){rownames(row_annot) = rownames(data)}
                          if(col_column){rownames(col_annot) = colnames(data)}
                          if(is.na(fontsize_row)){fontsize_row = cellheight}
                          if(is.na(fontsize_col)){fontsize_col = cellwidth}

                          OUTPUT %>% dirname() %>% dir.create_p()
                          
                          pheatmap::pheatmap(data,
                                             cluster_rows = cluster_rows, cluster_cols = cluster_cols , 
                                             cellwidth = cellwidth, cellheight = cellheight, color=heatmap_col,
                                             filename = OUTPUT, border_color = border_color ,na_col = "gray90",
                                             annotation_colors = annot_palette,
                                             annotation_row    = row_annot,
                                             annotation_col    = col_annot,
                                             gaps_row = gaps_row_num,
                                             gaps_col = gaps_col_num,
                                             display_numbers = display_numbers,
                                             fontsize_number = fontsize_number, number_color = "black",
                                             fontsize = fontsize , fontsize_row = fontsize_row, fontsize_col = fontsize_col, 
                                             cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                                             clustering_method = clustering_method,
                                             clustering_distance_rows = clustering_distance_rows,
                                             clustering_distance_cols = clustering_distance_cols,
                                             main = TITLE                                                      )
}

########################################################
################# Seurat assay heatmap   ###############
########################################################

assay_heatmap = function(Seurat        = ASAP,
                         assay         = "motif_zscore",
                         slot          = "data",
                         select_list   = NA,
                         meta_name     = "HTO_final_4",
                         target_list   = HTO_list_2,
                         upper_thresh  = 2,
                         lower_thresh  = -2,
                         sampling      = TRUE,
                         sample_num    = 10000,
                         TITLE         = "title",
                         OUTPUT        = "tmp.pdf",
                         save_raw_data = TRUE,
                         OUTPUT_data   = "tmp.txt",
                         plot_h =20, plot_w =60
                         ){
                        data      = GetAssayData(Seurat, assay = assay,  slot = slot )  %>% t() %>% as.data.frame()
                        data      = data[,!apply(data,2,anyNA)]
                        if(!is.na(upper_thresh)){data[data>upper_thresh] = upper_thresh}
                        if(!is.na(lower_thresh)){data[data<lower_thresh] = lower_thresh}
                        if(!is.na(select_list)[1]){data = data[,select_list]}

                        meta_tmp  = Seurat@meta.data[,meta_name,drop=F] %>%
                                     dplyr::rename(target=1)
                        
                        tmp_list = list()
                        for(tmp_num in 1:nrow(target_list)){
                            paste0(tmp_num," : ",target_list$HTO[tmp_num]) %>% print()
                            tmp_CB_num = sum(meta_tmp$target == target_list$HTO[tmp_num])
                            if(tmp_CB_num != 0){CellBarcode_tmp    = meta_tmp %>% dplyr::filter(target == target_list$HTO[tmp_num]) %>% rownames(.) }
                        
                            if(tmp_CB_num >1 ){tmp_list[[tmp_num]] = hclust_row_order(data[CellBarcode_tmp,],METHOD="ward.D2")}
                            if(tmp_CB_num ==1){tmp_list[[tmp_num]] = CellBarcode_tmp }
                        }
                        
                        name_list   = unlist(tmp_list)

                        if(sampling){
                          sampling_tmp = sample(1:nrow(data), min(sample_num,nrow(data)), replace = FALSE) %>% sort()
                          name_list_2 = hclust_col_order(data[sampling_tmp,],METHOD="ward.D2")
                          }else{
                          name_list_2 = hclust_col_order(data,METHOD="ward.D2")
                        }
                        
                        df = data[name_list,name_list_2] %>% t() %>% as.data.frame() %>%
                              rownames_to_column("rownames") %>% gather(.,key=colnames,value=score,-rownames)
                        if(save_raw_data){write.table_n_2(data[name_list,name_list_2],"CellBarcode",OUTPUT_data)}

                        df2 = left_join(data.frame(CellBarcode=df$colnames%>%unique(.)) , meta_tmp %>% rownames_to_column("CellBarcode"),by="CellBarcode") %>%
                               mutate(target = factor(target,levels=target_list$HTO) )
                        
                        tmp_target_num  = rep(NA,nrow(target_list))
                        for(i in 1:nrow(target_list)){tmp_target_num[i]=grep(paste0("^",target_list$HTO[i],"$"),df2$target)%>%max();names(tmp_target_num)[i]=target_list$HTO[i]}
                        
                        p1  = ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                            ggrastr::geom_tile_rast()+
                            scale_fill_gradientn(colours=PurpleAndYellow())+
                            plot_theme() +
                            scale_x_discrete(limits=unique(df$colnames)) +
                            scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
                            plot_theme(theme="void",legend=TRUE)+
                            theme(axis.text.x=element_blank(),
                                  strip.text=element_text(size=8),
                                  axis.line = element_line(colour = "black"),
                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                            labs(x="Cell Barcode",y="")
                        
                        for(i in 1:nrow(target_list)){
                          p1 = p1 +
                               geom_vline(xintercept=tmp_target_num[i]+0.5,col="black")
                        }
                        
                        p2 = ggplot(df2, aes(x=CellBarcode,y=1,fill=target))+
                              ggrastr::geom_tile_rast()+
                              scale_fill_manual(values = target_list$color)+
                              scale_y_continuous(expand=c(0,0)) +
                              scale_x_discrete(limits=unique(df$colnames)) +
                              theme(axis.title.x=element_blank(),axis.ticks=element_blank(),
                                    axis.text=element_blank(),
                                    axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
                                    plot.title = element_text(size=10,face="bold", family = "Helvetica"),
                                    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
                              labs(title=paste0(TITLE),y="target")
                        
                        p_sum = p2+p1+  plot_layout(ncol = 1, heights = c(1, 9))
                        
                        pdf_2(OUTPUT ,h=plot_h,w=plot_w)
                         plot(p_sum)
                        dev.off()
                      }

########################################################
#################  plot pie_chart        ###############
########################################################

plot_piechart = function(target_list,
                         category_list,
                         category_level,
                         category_color = NA,
                         TITLE = "TITLE",
                         OUTPUT = "scatter_pie/scattepie.pdf",
                         w_enlarge = 1,
                         h_enlarge = 1
                         ){
                           data = table(target_list,category_list)
                           data_2 = data %>% as.matrix.data.frame(.) %>% as.data.frame(.)
                           rownames(data_2) = rownames(data)
                           colnames(data_2) = colnames(data)
                           data_2$Size = rowSums(data_2)
                           ncol_tmp = sqrt(nrow(data_2) * 2 ) %>% ceiling(.)
                           nrow_tmp = (nrow(data_2)/ncol_tmp) %>% ceiling(.)
                           data_2[is.na(data_2)] = 0
                           data_2 = data_2 %>%
                                             mutate(x      = rep(1:ncol_tmp,nrow_tmp)      %>% head(nrow(data_2)),
                                                    y      = rep(nrow_tmp:1,each=ncol_tmp) %>% head(nrow(data_2)),
                                                    radius = log10(Size)/max(log10(Size)) * 0.4                            ) %>%
                                             rownames_to_column("Name")
                           
                           p = ggplot(data_2, aes(x = x , y = y )) +
                                scatterpie::geom_scatterpie(aes(x = x , y = y , r = radius),
                                                            data = data_2, cols = category_level, color = NA) +
                                theme_void() +
                                geom_point(size=6,shape=16,col="white")+
                                geom_text(aes(x = x, y = y + 0.5, label = Name),  size =3, ,fontface = "bold", family = "Helvetica", inherit.aes = FALSE)+
                                geom_text(aes(x = x, y = y , label = Size), size =2, ,fontface = "bold", family = "Helvetica", inherit.aes = FALSE)+
                                plot_theme(legend=TRUE,theme="void") +
                                labs(title=TITLE)+
                                theme(axis.text=element_blank(),axis.title=element_blank())
                           
                           if(!is.na(category_color[1])){p = p + scale_fill_manual(values = category_color) }
                           
                           paste0(OUTPUT) %>% pdf_2(.,h = nrow_tmp*h_enlarge , w = ncol_tmp*w_enlarge)
                            plot(p)
                           dev.off()
                          }

######################
sep_index = function(list,list_num=NA,sep=NA
                    ){
                      num_of_separate = length(list)
                      if(is.na(sep)){sep = ceiling( num_of_separate / list_num ) }
                      if(is.na(list_num)){list_num = ceiling( num_of_separate / sep ) }
                      tmp_list_num = num_of_separate- list_num * (sep - 1)
                      if(tmp_list_num!=list_num){sep_list = c(rep(1:tmp_list_num,each=sep),rep((tmp_list_num+1):list_num,each=(sep-1)))
                                           }else{sep_list = c(rep(1:list_num,each=sep))}
                      tmp_res = data.frame(target = list) %>%
                                 mutate(sep=sep_list)
                      tmp_res_2 = split(tmp_res,tmp_res$sep) %>%
                                   lapply(.,function(x){x$target})
                      names(tmp_res_2) = c()
                      tmp_res_2
                     }

######################
plot_colordot_rast = function(data,
                             color_high = "red",
                             color_low = "blue",
                             color_zero = "white",
                             colname="", rowname="",
                             title="title",
                             max_size=4,
                             size=10,
                             legend=TRUE,
                             theme="classic",
                             x_angle=90
                            ){data_2 = data %>%
                                        rownames_to_column("rowname") %>%
                                        pivot_longer(., col = -rowname, names_to = "colname", values_to = "value") %>%
                                        as.data.frame(.) %>%
                                        mutate(x=colname %>% factor(., levels=colnames(data))         %>% as.numeric(.),
                                               y=rowname %>% factor(., levels=rownames(data)%>%rev()) %>% as.numeric(.),
                                               size=abs(value))
                              ncol_tmp = ncol(data)
                              nrow_tmp = nrow(data)
                              
                              p1 = ggplot(data_2,aes(x=x,y=y,fill=value))
                              for(tmp_index in 1:(ncol_tmp-1)){
                                p1 = p1 + geom_vline(xintercept=tmp_index, col="grey", size=0.5, alpha=0.2, linetype="dashed") 
                              }
                              for(tmp_index in 1:(nrow_tmp-1)){
                                p1 = p1 + geom_hline(yintercept=tmp_index, col="grey", size=0.5, alpha=0.2, linetype="dashed") 
                              }
                              p1 = p1 +
                                    geom_point(aes(size=size), color="black", shape=21 )+
                                    scale_fill_gradient2(low=color_low, high=color_high, mid=color_zero, midpoint = 0)+
                                    scale_size_area(max_size=max_size, limits = c(0, 4) ) + 
                                    plot_theme(size=size, theme=theme, legend=legend, x_angle=x_angle) +
                                    scale_x_continuous(breaks=1:ncol_tmp,
                                                       labels=colnames(data),
                                                       expand=c(0,0),limits=c(0.4,ncol_tmp+0.6))+
                                    scale_y_continuous(breaks=1:nrow_tmp,
                                                       labels=rownames(data)%>%rev(.),
                                                       expand=c(0,0),limits=c(0.4,nrow_tmp+0.6))+
                                    labs(x=colname ,y=rowname, title=title )
                              p1
                             }

######################
HTODemux_f = function(object,assay = "CROP", positive.quantile = 0.999, max_cutoff_1st=20,
                     plot_hist=TRUE, hist_path=paste0("tmp.pdf"),
                     seed = 111, verbose = TRUE){
    ## https://github.com/satijalab/seurat/blob/HEAD/R/preprocessing.R
    ## if (!is.null(x = seed)) { set.seed(seed = seed) }
    #initial clustering
    ## assay <- assay %||% DefaultAssay(object = object)
    counts_tmp = GetAssayData(object = object, assay = assay, slot = 'counts' )
    data_tmp   = GetAssayData(object = object, assay = assay, slot = 'data' )
    discrete   = GetAssayData(object = object, assay = assay, slot = 'counts' )
    discrete[discrete > 0] <- 0
    
    single_num_list = list()
    for(tmp_cutoff in 1:max_cutoff_1st){
        cutoff <- rep(tmp_cutoff,nrow(data_tmp))
        names(cutoff) = rownames(x = data_tmp)
        
        # for each HTO, we will use the minimum cluster for fitting
        discrete_tmp = discrete
        for (iter in rownames(x = data_tmp)) {
          values <- counts_tmp[iter, colnames(object)]
          discrete_tmp[iter, names(x = which(x = values > cutoff[iter] ))] <- 1
        }
        
        single_num_list[[as.character(tmp_cutoff)]] = colSums(discrete_tmp) %>% table() %>% as.data.frame() %>%
                                                       dplyr::rename(target_num=1) %>%
                                                       dplyr::filter(target_num==1) %>%
                                                       pull(Freq)
        paste0("cutoff 1stStep : Cutoff ",tmp_cutoff," has finished.") %>% print()
    }
    
    single_num_final = names(single_num_list)[as.numeric(single_num_list) %>% which.max() %>% min()] %>% as.numeric_f()
    paste0("cutoff 1stStep : Cutoff : ", single_num_final," was best.") %>% print()
    if(single_num_final>=max_cutoff_1st){print("Error! Please Set the Cutoff score higher.")}
    
    for (iter in rownames(x = data_tmp)) {
      values <- counts_tmp[iter, colnames(object)]
      discrete[iter, names(x = which(x = values > single_num_final ))] <- 1
    }
    
    colSums(discrete) %>% table() %>% print()
    ##     0     1     2     3     4     5     6     7     8     9    10    11    12 
    ##  2045 10799  4361  1548   602   275   111    64    43    26    10    11     5 
    ##    13    14    15    22 
    ##     4     3     2     1 
    
    ## Definition of function
    ## https://github.com/satijalab/seurat/blob/master/R/utilities.R
    MaxN <- function(x, N = 2){
      len <- length(x)
      if (N > len) {
        warning('N greater than length(x).  Setting N=length(x)')
        N <- length(x)
      }
      sort(x, partial = len - N + 1)[len - N + 1]
    }
    
    # now assign cells to HTO based on discretized values
    npositive <- colSums(x = discrete)
    classification.global <- npositive
    classification.global[npositive == 0] <- "Negative"
    classification.global[npositive == 1] <- "Singlet"
    classification.global[npositive > 1] <- "Doublet"
    donor.id = rownames(x = data_tmp)
    hash.max <- apply(X = data_tmp, MARGIN = 2, FUN = max)
    hash.maxID <- apply(X = data_tmp, MARGIN = 2, FUN = which.max)
    hash.second <- apply(X = data_tmp, MARGIN = 2, FUN = MaxN, N = 2)
    hash.maxID <- as.character(x = donor.id[sapply(
      X = 1:ncol(x = data_tmp),
      FUN = function(x) {
        return(which(x = data_tmp[, x] == hash.max[x])[1])
      }
    )])
    hash.secondID <- as.character(x = donor.id[sapply(
      X = 1:ncol(x = data_tmp),
      FUN = function(x) {
        return(which(x = data_tmp[, x] == hash.second[x])[1])
      }
    )])
    hash.margin <- hash.max - hash.second
    doublet_id <- sapply(
      X = 1:length(x = hash.maxID),
      FUN = function(x) {
        return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
      }
    )
    
    # doublet_names <- names(x = table(doublet_id))[-1] # Not used
    classification <- classification.global
    classification[classification.global == "Negative"] <- "Negative"
    classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
    classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
    classification.metadata <- data.frame(
      hash.maxID,
      hash.secondID,
      hash.margin,
      classification,
      classification.global
    )
    
    colnames(x = classification.metadata) <- paste(
      assay,"1stStep",
      c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
      sep = '_'
    )
    object <- AddMetaData(object = object, metadata = classification.metadata)
    
    doublets <- rownames(x = object[[]])[which(object[[paste(assay,"1stStep", "classification.global",sep = '_')]] == "Doublet")]
    Idents(object) <- paste(assay,"1stStep", 'classification',sep = '_')
    Idents(object = object, cells = doublets) <- 'Doublet'
    print("cutoff 1stStep finished.")
    
    ################################
    
    ## 2nd step 
    ## For certain conditions, CBs marked as the other condition in the first step are set as the negative control.  (instead of k-means/clara cluster)
    discrete <- GetAssayData(object = object, assay = assay)
    discrete[discrete > 0] <- 0
    
    # for each HTO, we will use the minimum cluster for fitting
    cut_off_list = data.frame(condition=rownames(x = data_tmp), first_step=single_num_final, second_step=NA) %>%
                    column_to_rownames("condition")
    for (iter in rownames(x = data_tmp)) {
      tmp_CB = names(Idents(object))[Idents(object) %in% setdiff(rownames(x = data_tmp),iter)]
      values_use <- counts_tmp[iter, tmp_CB]
      values     <- counts_tmp[iter, ]
    
      fit <- suppressWarnings(expr = fitdistrplus::fitdist(data = values_use, distr = "nbinom"))
      cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
      discrete[iter, names(x = which(x = values > cutoff))] <- 1
      cut_off_list[iter,"second_step"]=cutoff
      if (verbose) {
        message(paste0("Cutoff 2nd step for ", iter, " : ", cutoff, " reads"))
      }
    }
    
    ## plot histogram
    if(plot_hist){
        p = list()
        for (iter in rownames(x = data_tmp)) {
          values <- counts_tmp[iter, colnames(object)]
        
          df = data.frame(CB= names(values),value=values) %>%
                mutate(value2 = ifelse(value>30,30,value))
          p[[iter]] = ggplot(df, aes(x=value2))+
                       ## geom_density(alpha = .3, col = "black",fill="green") +
                       geom_histogram(binwidth = 1,alpha=0.5)+
                       scale_y_log10(expand=c(0,0),breaks=10^(0:4),limits = c(1,20000) ,labels=scales::trans_format("log10",scales::math_format(10^.x)))+
                       geom_vline(xintercept=cut_off_list[iter,"first_step"]+0.5,col="blue",size=0.5)+
                       geom_vline(xintercept=cut_off_list[iter,"second_step"]+0.5,col="red",size=0.5)+
                       plot_theme()+
                       labs(title=paste0(iter," : ",sum(df$value2 > cut_off_list[iter,"second_step"]+0.5 ),"Cells , 1st:",cut_off_list[iter,"first_step"],"reads, 2nd:",cut_off_list[iter,"second_step"],"reads"), x="counts")
        }
        
        ncol = 20
        p_sum = patchwork::wrap_plots(p, ncol)+
                 plot_annotation("Counts gRNA in each condition")
        nrow = length(p)/ncol %>% ceiling(.)
        
        hist_path %>% pdf_2(.,h=5*nrow,w=5*ncol)
         plot(p_sum)
        dev.off()
    }
    
    # now assign cells to HTO based on discretized values
    npositive <- colSums(x = discrete)
    classification.global <- npositive
    classification.global[npositive == 0] <- "Negative"
    classification.global[npositive == 1] <- "Singlet"
    classification.global[npositive > 1] <- "Doublet"
    donor.id = rownames(x = data_tmp)
    hash.max <- apply(X = data_tmp, MARGIN = 2, FUN = max)
    hash.maxID <- apply(X = data_tmp, MARGIN = 2, FUN = which.max)
    hash.second <- apply(X = data_tmp, MARGIN = 2, FUN = MaxN, N = 2)
    hash.maxID <- as.character(x = donor.id[sapply(
      X = 1:ncol(x = data_tmp),
      FUN = function(x) {
        return(which(x = data_tmp[, x] == hash.max[x])[1])
      }
    )])
    hash.secondID <- as.character(x = donor.id[sapply(
      X = 1:ncol(x = data_tmp),
      FUN = function(x) {
        return(which(x = data_tmp[, x] == hash.second[x])[1])
      }
    )])
    hash.margin <- hash.max - hash.second
    doublet_id <- sapply(
      X = 1:length(x = hash.maxID),
      FUN = function(x) {
        return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), collapse = "_"))
      }
    )
    
    # doublet_names <- names(x = table(doublet_id))[-1] # Not used
    classification <- classification.global
    classification[classification.global == "Negative"] <- "Negative"
    classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
    classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
    classification.metadata <- data.frame(
      hash.maxID,
      hash.secondID,
      hash.margin,
      classification,
      classification.global
    )
    
    colnames(x = classification.metadata) <- paste(
      assay,"2ndStep",
      c('maxID', 'secondID', 'margin', 'classification', 'classification.global'),
      sep = '_'
    )
    object <- AddMetaData(object = object, metadata = classification.metadata)
    
    doublets <- rownames(x = object[[]])[which(object[[paste(assay,"2ndStep", "classification.global",sep = '_')]] == "Doublet")]
    Idents(object) <- paste(assay,"2ndStep", 'classification',sep = '_')
    Idents(object = object, cells = doublets) <- 'Doublet'
    print("cutoff 2ndStep finished.")
    return(object)
  }

########################################################
#################         Signac         ###############
########################################################

library(Signac)

FindRegion2 <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

comet = c("#E6E7E8","#3A97FF","#8816A7","black")
greenBlue = c('#e0f3db','#ccebc5','#a8ddb5','#4eb3d3','#2b8cbe','#0868ac','#084081')

LinkPlot2 <- function(
  object,
  region,
  assay = NULL,
  min.cutoff = 0,
  sep = c("-", "-"),
  extend.upstream = 0,
  extend.downstream = 0,
  scale.linewidth = FALSE
) {
  region <- FindRegion2(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  chromosome <- seqnames(x = region)

  # extract link information
  links <- Links(object = object)

  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }

  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)

  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  link.df <- link.df[abs(x = link.df$score) > min.cutoff, ]

  # remove links outside region
  link.df <- link.df[link.df$start >= start(x = region) & link.df$end <= end(x = region), ]

  # plot
  if (nrow(x = link.df) > 0) {
    if (!requireNamespace(package = "ggforce", quietly = TRUE)) {
      warning("Please install ggforce to enable LinkPlot plotting: ",
              "install.packages('ggforce')")
      p <- ggplot(data = link.df)
    } else {
      # convert to format for geom_bezier
      link.df$group <- seq_len(length.out = nrow(x = link.df))
      df <- data.frame(
        x = c(link.df$start,
              (link.df$start + link.df$end) / 2,
              link.df$end),
        y = c(rep(x = 0, nrow(x = link.df)),
              rep(x = -1, nrow(x = link.df)),
              rep(x = 0, nrow(x = link.df))),
        group = rep(x = link.df$group, 3),
        score = rep(link.df$score, 3)
      )
      min.color <- min(0, min(df$score))
      if (scale.linewidth) {
        p <- ggplot(data = df) +
          ggforce::geom_bezier(
            mapping = aes_string(x = "x", y = "y", group = "group", color = "score", linewidth = "score")
          )
      } else {
        p <- ggplot(data = df) +
          ggforce::geom_bezier(
            mapping = aes_string(x = "x", y = "y", group = "group", color = "score")
          )
      }
      p <- p +
        geom_hline(yintercept = 0, color = 'grey') +
    #    scale_color_gradient2(low = "#f1eef6", mid = "#d4b9da", high = "#ce1256",
    #                          limits = c(min.color, max(df$score)),
    #                          n.breaks = 5)
        scale_colour_gradientn(colours = greenBlue, limits = c(0, max(df$score)))
    }
  } else {
    p <- ggplot(data = link.df)
  }
  p <- p +
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    ylab("Links") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
  return(p)
}

####################
Extract_bed_path <- function(TF,lineage){
    dir <- file.path("Reproducibility","Results","motifmatchr", lineage)
    hits <- list.files(dir, pattern = paste0("^", TF, ".*_peak\\.bed$"),
                       full.names = TRUE)
    if (length(hits) == 0) stop("No BED found for ", TF, " in ", dir)
    if (length(hits) > 1) {
      hits <- hits[which.max(file.info(hits)$mtime)]
    }
    bed_path <- hits
    bed_path
}

####################
Plot_genetrack <- function(Obj,group.by,lineage,group_list,gene,region_list,motif_list,cutoff,cols,path){
    region = region_list[[gene]]
    motifs = motif_list[[gene]]

    # Extract the chromosome, start, and end positions
    gr = StringToGRanges(region)
    chromosome <- as.character(seqnames(gr))
    start_position <- start(gr)
    end_position <- end(gr)

    # load links
    links = file.path("Reproducibility","Results","Cicero", lineage, paste0("UC_DOGMA_links_",lineage,"_",chromosome, ".rds")) %>% readRDS()
    cA = links %>% as.data.frame()
    cA_tmp = dplyr::filter(cA, score > cutoff)

    DefaultAssay(Obj) = 'ATAC'
    Links(Obj) <- GRanges(ranges=IRanges(start=cA_tmp$start, end=cA_tmp$end), 
                            seqnames=cA_tmp$seqnames,
                            score=cA_tmp$score)

    # Define the tmp_region as a GRanges object
    bed_list = list()
    for(tmp_motif in motifs){
        bed_list[[tmp_motif]] <- Extract_bed_path(TF=tmp_motif,lineage=lineage) %>%
                 fread() %>% as.data.frame() %>% dplyr::mutate(center = (V2+V3)/2) %>%
                 dplyr::mutate(range = paste0(V1,"-",center-400,"-",center+400)) %>% pull("range")
    }

    bed_total = unlist(bed_list)
    ranges.show <- subsetByOverlaps(StringToGRanges(bed_total), gr)
    ranges.show$color <- "orange"

    cov_plot <- CoveragePlot(
        object = Obj,
        assay = 'ATAC',
        region = region,
        window = 500,
        peaks = FALSE,
        links = FALSE,
        annotation = FALSE,
        region.highlight = ranges.show,
        idents = group_list,
        group.by=group.by
    )

    gene_plot <- AnnotationPlot(
      object = Obj,
      region = region
    )

    link_plot <- LinkPlot2(
        object = Obj,
        region = region,
        min.cutoff = cutoff
    )

    motif_plots <- imap(bed_list, function(ranges, TF) {
      if (is.null(ranges)) return(NULL)
      PeakPlot(
        object = Obj,
        region = region,
        peaks  = StringToGRanges(ranges),
        color  = "brown3"
      ) + ylab(TF)
    }) |> compact()

    expr_plot= ExpressionPlot(
        object = Obj,
        features = gene,
        group.by = group.by,
        assay = "RNA",
        idents = group_list,
        )
    
    ## 3) Combine tracks
    plotlist <- c(list(cov_plot, gene_plot, link_plot), motif_plots)
    heights  <- c(8, 3, 2, rep(0.5, length(motif_plots)))
    
    p <- CombineTracks(
      plotlist        = plotlist,
      expression.plot = expr_plot,
      heights         = heights,
      widths          = c(10, 1)
    )
    p = p & scale_fill_manual(values = cols)
    paste0(path) %>% pdf_3(., w=12,h=6)
     print(p)
    dev.off()
}

####################
Plot_genetrack2 <- function(Obj,lineage,celltype_list,gene,region_list,motif_list,cutoff,cols,path){
    region = region_list[[gene]]
    motifs = motif_list[[gene]]

    # Extract the chromosome, start, and end positions
    gr = StringToGRanges(region)
    chromosome <- as.character(seqnames(gr))
    start_position <- start(gr)
    end_position <- end(gr)

    # load links
    links = file.path("Reproducibility","Results","Cicero", lineage, paste0("UC_DOGMA_links_",lineage,"_",chromosome, ".rds")) %>% readRDS()
    cA = links %>% as.data.frame()
    cA_tmp = dplyr::filter(cA, score > cutoff)

    DefaultAssay(Obj) = 'ATAC'
    Links(Obj) <- GRanges(ranges=IRanges(start=cA_tmp$start, end=cA_tmp$end), 
                            seqnames=cA_tmp$seqnames,
                            score=cA_tmp$score)

    # Define the tmp_region as a GRanges object
    bed_list = list()
    for(tmp_motif in motifs){
        bed_list[[tmp_motif]] <- Extract_bed_path(TF=tmp_motif,lineage=lineage) %>%
                 fread() %>% as.data.frame() %>% dplyr::mutate(center = (V2+V3)/2) %>%
                 dplyr::mutate(range = paste0(V1,"-",center-400,"-",center+400)) %>% pull("range")
    }

    bed_total = unlist(bed_list)
    ranges.show <- subsetByOverlaps(StringToGRanges(bed_total), gr)
    ranges.show$color <- "orange"

    cov_plot <- CoveragePlot(
        object = Obj,
        assay = 'ATAC',
        region = region,
        window = 500,
        peaks = FALSE,
        links = FALSE,
        annotation = FALSE,
        region.highlight = ranges.show,
        idents = celltype_list,
        group.by='celltype'
    )

    gene_plot <- AnnotationPlot(
        object = Obj,
        region = region,
        mode = "transcript"
    )

    link_plot <- LinkPlot2(
        object = Obj,
        region = region,
        min.cutoff = cutoff
    )

    motif_plots <- imap(bed_list, function(ranges, TF) {
      if (is.null(ranges)) return(NULL)
      PeakPlot(
        object = Obj,
        region = region,
        peaks  = StringToGRanges(ranges),
        color  = "brown3"
      ) + ylab(TF)
    }) |> compact()

    expr_plot= ExpressionPlot(
        object = Obj,
        features = gene,
        group.by = 'celltype',
        assay = "RNA",
        idents = celltype_list,
        )
    
    ## 3) Combine tracks
    plotlist <- c(list(cov_plot, gene_plot, link_plot), motif_plots)
    heights  <- c(8, 3, 2, rep(0.5, length(motif_plots)))
    
    p <- CombineTracks(
      plotlist        = plotlist,
      expression.plot = expr_plot,
      heights         = heights,
      widths          = c(10, 1)
    )
    p = p & scale_fill_manual(values = cols)
    paste0(path) %>% pdf_3(., w=12,h=6)
     print(p)
    dev.off()
}

####################
Plot_genetrack_nhood <- function(Obj,lineage,group_list,gene,region_list,motif_list,cutoff,cols,path){
    region = region_list[[gene]]
    motifs = motif_list[[gene]]

    # Extract the chromosome, start, and end positions
    gr = StringToGRanges(region)
    chromosome <- as.character(seqnames(gr))
    start_position <- start(gr)
    end_position <- end(gr)

    # load links
    links = file.path("Reproducibility","Results","Cicero", lineage, paste0("UC_DOGMA_links_",lineage,"_",chromosome, ".rds")) %>% readRDS()
    cA = links %>% as.data.frame()
    cA_tmp = dplyr::filter(cA, score > cutoff)

    DefaultAssay(Obj) = 'ATAC'
    Links(Obj) <- GRanges(ranges=IRanges(start=cA_tmp$start, end=cA_tmp$end), 
                            seqnames=cA_tmp$seqnames,
                            score=cA_tmp$score)

    # Define the tmp_region as a GRanges object
    bed_list = list()
    for(tmp_motif in motifs){
        bed_list[[tmp_motif]] <- Extract_bed_path(TF=tmp_motif,lineage=lineage) %>%
                 fread() %>% as.data.frame() %>% dplyr::mutate(center = (V2+V3)/2) %>%
                 dplyr::mutate(range = paste0(V1,"-",center-400,"-",center+400)) %>% pull("range")
    }

    bed_total = unlist(bed_list)
    ranges.show <- subsetByOverlaps(StringToGRanges(bed_total), gr)
    ranges.show$color <- "orange"

    cov_plot <- CoveragePlot(
        object = Obj,
        assay = 'ATAC',
        region = region,
        window = 500,
        peaks = FALSE,
        links = FALSE,
        annotation = FALSE,
        region.highlight = ranges.show,
        idents = group_list,
        group.by='nhood_groups'
    )

    gene_plot <- AnnotationPlot(
      object = Obj,
      region = region
    )

    link_plot <- LinkPlot2(
        object = Obj,
        region = region,
        min.cutoff = cutoff
    )

    motif_plots <- imap(bed_list, function(ranges, TF) {
      if (is.null(ranges)) return(NULL)
      PeakPlot(
        object = Obj,
        region = region,
        peaks  = StringToGRanges(ranges),
        color  = "brown3"
      ) + ylab(TF)
    }) |> compact()

    expr_plot= ExpressionPlot(
        object = Obj,
        features = gene,
        group.by = 'nhood_groups',
        assay = "RNA",
        idents = group_list,
        )
    
    ## 3) Combine tracks
    plotlist <- c(list(cov_plot, gene_plot, link_plot), motif_plots)
    heights  <- c(8, 3, 2, rep(0.5, length(motif_plots)))
    
    p <- CombineTracks(
      plotlist        = plotlist,
      expression.plot = expr_plot,
      heights         = heights,
      widths          = c(10, 1)
    )
    p = p & scale_fill_manual(values = cols)
    paste0(path) %>% pdf_3(., w=12,h=6)
     print(p)
    dev.off()
}

####################
Plot_genetrack_BCG <- function(Obj,lineage,group_list,gene,region_list,motif_list,cutoff,cols,path){
    region = region_list[[gene]]
    motifs = motif_list[[gene]]

    # Extract the chromosome, start, and end positions
    gr = StringToGRanges(region)
    chromosome <- as.character(seqnames(gr))
    start_position <- start(gr)
    end_position <- end(gr)

    # load links
    links = file.path("Reproducibility","Results","Cicero", lineage, paste0("UC_DOGMA_links_",lineage,"_",chromosome, ".rds")) %>% readRDS()
    cA = links %>% as.data.frame()
    cA_tmp = dplyr::filter(cA, score > cutoff)

    DefaultAssay(Obj) = 'ATAC'
    Links(Obj) <- GRanges(ranges=IRanges(start=cA_tmp$start, end=cA_tmp$end), 
                            seqnames=cA_tmp$seqnames,
                            score=cA_tmp$score)

    # Define the tmp_region as a GRanges object
    bed_list = list()
    for(tmp_motif in motifs){
        bed_list[[tmp_motif]] <- Extract_bed_path(TF=tmp_motif,lineage="BCG") %>%
                 fread() %>% as.data.frame() %>% dplyr::mutate(center = (V2+V3)/2) %>%
                 dplyr::mutate(range = paste0(V1,"-",center-400,"-",center+400)) %>% pull("range")
    }

    bed_total = unlist(bed_list)
    ranges.show <- subsetByOverlaps(StringToGRanges(bed_total), gr)
    ranges.show$color <- "orange"

    cov_plot <- CoveragePlot(
        object = Obj,
        assay = 'ATAC',
        region = region,
        window = 500,
        peaks = FALSE,
        links = FALSE,
        annotation = FALSE,
        region.highlight = ranges.show,
        idents = group_list,
        group.by='prepos'
    )

    gene_plot <- AnnotationPlot(
      object = Obj,
      region = region
    )

    link_plot <- LinkPlot2(
        object = Obj,
        region = region,
        min.cutoff = cutoff
    )

    motif_plots <- imap(bed_list, function(ranges, TF) {
      if (is.null(ranges)) return(NULL)
      PeakPlot(
        object = Obj,
        region = region,
        peaks  = StringToGRanges(ranges),
        color  = "brown3"
      ) + ylab(TF)
    }) |> compact()

    expr_plot= ExpressionPlot(
        object = Obj,
        features = gene,
        group.by = 'prepos',
        assay = "RNA",
        idents = group_list,
        )
    
    ## 3) Combine tracks
    plotlist <- c(list(cov_plot, gene_plot, link_plot), motif_plots)
    heights  <- c(8, 3, 2, rep(0.5, length(motif_plots)))
    
    p <- CombineTracks(
      plotlist        = plotlist,
      expression.plot = expr_plot,
      heights         = heights,
      widths          = c(10, 1)
    )
    p = p & scale_fill_manual(values = cols)
    paste0(path) %>% pdf_3(., w=12,h=6)
     print(p)
    dev.off()
}

####################
vlnplot_w_box <- function(Seurat, assay = 'RNA', slot = 'data', features, cols, file_name){
  plot_case <- function(signature, y_max = NULL){
    VlnPlot(Seurat, 
            features = signature,  # pass a single feature
            slot = slot,
            pt.size = 0, 
            group.by = "celltype", 
            cols = cols) + NoLegend() +
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)
  }

  plot_list <- list()
  
    for (tmp_feature in features) {
     plot_list[[tmp_feature]] <- plot_case(tmp_feature)
    }

   pdf(file_name, width = 4, height = 5)
     print(plot_list)
   dev.off()
 }
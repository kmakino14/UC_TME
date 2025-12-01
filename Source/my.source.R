suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(compositions))
options(stringsAsFactors=F)

take_factor <- function(list,order,sep){sapply(strsplit(list,sep),function(x){if(length(order)==1){x[order]}else{paste(x[order],collapse=sep)}})}

as.numeric_f<-function(number_list){as.numeric(as.character(number_list))}

last_col_names<-function(dataframe){colnames(dataframe)[ncol(dataframe)]}
last_row_names<-function(dataframe){rownames(dataframe)[nrow(dataframe)]}

take_same_column<-function(dataframe,ref_dataframe){
  dataframe[,intersect(colnames(dataframe),colnames(ref_dataframe))]}

correlation_p<-function(x,y){x<-as.numeric(x)
y<-as.numeric(y)
if(max(x)==min(x)){tmp<-c(rep("NA",4))}else if(max(y)==min(y)){tmp<-c(rep("NA",4))} else{
  c_result<-cor.test(x,y,method="pearson")
  tmp<-as.numeric(c(c_result$p.value,c_result$estimate,c_result$conf.int[1],c_result$conf.int[2]))}
names(tmp)<-c("cor_p_value","cor","cor_D95","cor_U95")
tmp}

correlation_s<-function(x,y){x<-as.numeric(x)
y<-as.numeric(y)
if(max(x)==min(x)){tmp<-c(rep("NA",2))}else if(max(y)==min(y)){tmp<-c(rep("NA",2))} else{
  c_result<-cor.test(x,y,method="spearman")
  tmp<-as.numeric(c(c_result$p.value,c_result$estimate))}
names(tmp)<-c("p_value_s","rs_s")
tmp}

regression<-function(x,y){x<-as.numeric(x)
y<-as.numeric(y)
if(max(x)==min(x)){tmp<-c(rep("NA",3))}else if(max(y)==min(y)){tmp<-c(rep("NA",3))}else{
  res<-lm(y~x)
  tmp<-as.numeric(cbind(summary(res)$coefficients[2,4],res$coefficients[2],res$coefficients[1]))}
names(tmp)<-c("reg_p_value","slope","Intercept")
tmp}

write.table_TT   <-function(dataframe,filename){write.table(dataframe,filename,row.names=T, col.names=T, sep="\t", append=F, quote=F)}
write.table_TT_2 <-function(dataframe,filename){dirname(filename)%>%dir.create_p();write.table(dataframe,filename,row.names=T, col.names=T, sep="\t", append=F, quote=F)}
write.table_TF   <-function(dataframe,filename){write.table(dataframe,filename,row.names=T, col.names=F, sep="\t", append=F, quote=F)}
write.table_TF_2 <-function(dataframe,filename){dirname(filename)%>%dir.create_p();write.table(dataframe,filename,row.names=T, col.names=F, sep="\t", append=F, quote=F)}
write.table_FT   <-function(dataframe,filename){write.table(dataframe,filename,row.names=F, col.names=T, sep="\t", append=F, quote=F)}
write.table_FT_2 <-function(dataframe,filename){dirname(filename)%>%dir.create_p();write.table(dataframe,filename,row.names=F, col.names=T, sep="\t", append=F, quote=F)}
write.table_FF   <-function(dataframe,filename){write.table(dataframe,filename,row.names=F, col.names=F, sep="\t", append=F, quote=F)}
write.table_FF_2 <-function(dataframe,filename){dirname(filename)%>%dir.create_p();write.table(dataframe,filename,row.names=F, col.names=F, sep="\t", append=F, quote=F)}
write.table_n    <-function(dataframe,name_of_row,filename){dataframe = dataframe %>% rownames_to_column(name_of_row);write.table(dataframe,filename,row.names=F, col.names=T, sep="\t", append=F, quote=F)}
write.table_n_2  <-function(dataframe,name_of_row,filename){dirname(filename)%>%dir.create_p();dataframe = dataframe %>% rownames_to_column(name_of_row);write.table(dataframe,filename,row.names=F, col.names=T, sep="\t", append=F, quote=F)}

read.table_TT<-function(filename){read.table(filename,header=T,row.names=1,stringsAsFactors=F)}
read.table_TF<-function(filename){read.table(filename,header=F,row.names=1,stringsAsFactors=F)}
read.table_FT<-function(filename){read.table(filename,header=T,row.names=NULL,stringsAsFactors=F)}
read.table_FF<-function(filename){read.table(filename,header=F,row.names=NULL,stringsAsFactors=F)}

read.table_TT_sepT<-function(filename){read.table(filename,header=T,row.names=1,stringsAsFactors=F,sep="\t")}
read.table_TF_sepT<-function(filename){read.table(filename,header=F,row.names=1,stringsAsFactors=F,sep="\t")}
read.table_FT_sepT<-function(filename){read.table(filename,header=T,row.names=NULL,stringsAsFactors=F,sep="\t")}
read.table_FF_sepT<-function(filename){read.table(filename,header=F,row.names=NULL,stringsAsFactors=F,sep="\t")}

read.csv_TT = function(filename){read.table(filename,header=T,row.names=1,stringsAsFactors=F,sep=",")}
read.csv_TF = function(filename){read.table(filename,header=F,row.names=1,stringsAsFactors=F,sep=",")}
read.csv_FT = function(filename){read.table(filename,header=T,row.names=NULL,stringsAsFactors=F,sep=",")}
read.csv_FF = function(filename){read.table(filename,header=F,row.names=NULL,stringsAsFactors=F,sep=",")}

pdf_2 = function(filename){dirname(filename)%>%dir.create_p();pdf(filename)}
pdf_3 = function(filename,h=h,w=w){dirname(filename)%>%dir.create_p();pdf(filename,h=h,w=w)}
png_2 = function(filename){dirname(filename)%>%dir.create_p();png(filename)}
png_3 = function(filename,h=h,w=w){dirname(filename)%>%dir.create_p();png(filename,h=h,w=w)}

pick_up_file<-function(dir,letter){dir(dir)[grep(letter,dir(dir))]}

exchange_col<-function(dataframe,ncol1,ncol2){tmp1<-dataframe[,ncol1]
tmp2<-dataframe[,ncol2]
dataframe[,ncol1]<-tmp2
dataframe[,ncol2]<-tmp1
dataframe}

exchange_row<-function(dataframe,nrow1,nrow2){tmp1<-dataframe[nrow1,]
tmp2<-dataframe[nrow2,]
dataframe[nrow1,]<-tmp2
dataframe[nrow2,]<-tmp1
dataframe}

is.integer <-function(x){x%%1==0}
is.integer2<-function(x){ifelse(is.integer(x),x,NA)}

dir.create_p<-function(dirname){tmp<-strsplit(dirname,"/")
tmp<-tmp[[1]]
tmp2<-length(tmp)
name<-tmp[1]
for (i in 1:tmp2){if(file.exists(name)==1){cat("")} else(dir.create(name))
  name<-paste(name,tmp[i+1],sep="/")}   }

today<-function(){tmp<-strsplit(as.character(Sys.Date()),"-")[[1]]
                  paste0(substr(tmp[1],3,4),tmp[2],tmp[3])}

paste_u<-function(...){paste (..., sep = "_", collapse = NULL)}
paste_c<-function(...){paste (..., sep = "-", collapse = NULL)}
paste_s<-function(...){paste (..., sep = "/", collapse = NULL)}
paste_tmp = function(list){paste(list,collapse="__")}

change_col_names<-function(dataframe,col_name_old,col_name_new){colnames(dataframe)[match(col_name_old,colnames(dataframe))]<-col_name_new;dataframe}
before_write_table<-function(dataframe,name_of_row){dataframe<-cbind(rownames(dataframe),dataframe);colnames(dataframe)[1]<-name_of_row;dataframe}

fread_n   = function(filename){data<-fread(filename);data<-as.data.frame(data);rownames(data)<-data[,1];data<-data[,-1,drop=F];data}
fread_FT  = function(filename){data<-fread(filename);data<-as.data.frame(data);data}
fread_FF  = function(filename){data<-fread(filename);data<-as.data.frame(data);data=rbind(colnames(data),data);colnames(data)=paste0("V",1:ncol(data));data}
## fread_TT = function(filename){data<-fread(filename);data<-as.data.frame(data);rownames(data)<-data[,1];data<-data[,-1,drop=F];data}

pick_up_col<-function(dataframe,target){dataframe[,grep(target,colnames(dataframe))]}
pick_up_row<-function(dataframe,target){dataframe[grep(target,colnames(dataframe)),]}

merge_list<-function(list){paste0("^",paste0(list,collapse="$|^"),"$")}

cpm_n<-function(vector){vector/sum(vector)*10^6}

max_normalized<-function(each_row){each_row/max(each_row)}

list.files_n<-function(dir,pattern){list.files(dir,pattern=pattern,recursive=T,full.names=T)}

ggColorHue <- function(n, l=65) {hues <- seq(15, 375, length=n+1);hcl(h=hues, l=l, c=100)[1:n]}

q <- function (save = "no", status = 0, runLast = TRUE){.Internal(quit(save, status, runLast))}

p_convert      = function(x){if(is.na(x)){NA}else{as.character(symnum(x, cut=c(1,0.05,0.01,0.001,0.0001,0.00001,0.0000005,0), sym=c("******","*****","****","***","**","*","")))}}
thresh_convert = function(x,Up_limit=1000,thresh=10,Down_limit=0,if_up=1,if_down=0){if(is.na(x)){y=x}else{y=symnum(x, cut=c(Up_limit,thresh,Down_limit), sym=c(if_down,if_up));as.numeric_f(y)}}

log10P_convert    = function(data.frame){apply(data.frame,2,function(x){x[x>0.05]=1;y=-log10(x);y})}
log10P_convert_2  = function(data.frame){apply(data.frame,2,function(x){x[is.na(x)]=1;y=-log10(x);y})}
log10P_convert_v  = function(vector){vector[vector>0.05]=1;y=-log10(vector);y}
# log10P_convert_NA = function(data.frame){apply(data.frame,2,function(x){x[is.na(x)]=1;x[x>0.05]=1;y=-log10(x);y})}

ScientificNotation <- function(l) {l <- format(l, scientific = TRUE)
     l <- gsub("^(.*)e", "'\\1'e", l)
     l <- gsub("e\\+", "%*%10^", l)
     l[1] <- "0"
     parse(text = l)  }


SN_convert_2=function(num,digit,sep){
  x  = num
  digit=digit
  x1 = format(x,digits=digit,scientific = TRUE)
  x2 = take_factor(x1,1,"e") %>% as.numeric_f() %>% formatC(., digits = digit-1, format = "f")
  x3 = take_factor(x1,2,"e") %>% as.numeric_f()
  x4 = if(x3>=0){formatC(x3, width=2,flag="0")}else{formatC(x3, width=3,flag="0")}
  x5 = paste0(x2,sep,x4)
  x5
}

SN_convert = function(num,digit=2,sep="xE"){if(abs(num)<0.01){SN_convert_2(num,digit=digit,sep=sep)}else{round(num,digit=2)}}

formatC_2 = function(x,w=2){formatC(x,width=w,flag="0")}


PTscore <- function(obj, txs,
                     seqlev=intersect(seqlevels(obj), seqlevels(txs)),
                     upstream=2000, downstream=500){
  stopifnot(is(obj, "GAlignments"))
  stopifnot(is(txs, "GRanges"))
  obj <- as(obj, "GRanges")
  mcols(obj) <- NULL
  obj <- promoters(obj, upstream = 0, downstream = 1)
  cvg <- coverage(obj)
  cvg <- cvg[sapply(cvg, mean)>0]
  cvg <- cvg[names(cvg) %in% seqlev]
  seqlev <- seqlev[seqlev %in% names(cvg)]
  cvg <- cvg[seqlev]
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  pro <- promoters(txs, upstream = upstream, downstream = downstream)
  body <- shift(pro, shift = upstream + downstream)
  body[strand(pro)=="-"] <- shift(body[strand(pro)=="-"],
                                  shift = -1 * (upstream + downstream))
  pro$source <- "promoter"
  body$source <- "transcript"
  pro$oid <- seq_along(pro)
  body$oid <- seq_along(body)

  sel.gr <- c(pro, body)
  sel.gr <- sel.gr[order(sel.gr$oid)]

  sel.gr <- split(sel.gr, seqnames(sel.gr))
  seqlev <- seqlev[seqlev %in% names(sel.gr)]
  sel.gr <- sel.gr[seqlev]
  cvg <- cvg[seqlev]
  vws <- Views(cvg, sel.gr)
  vms <- viewMeans(vws)
  sel.gr <- unlist(sel.gr)
  sel.gr$score <- unlist(vms)

  pro <- sel.gr[sel.gr$source %in% "promoter"]
  body <- sel.gr[sel.gr$source %in% "transcript"]
  stopifnot(identical(pro$oid, body$oid))
  sel <- txs
  sel$promoterPart <- ranges(pro)
  sel$transcriptPart <- ranges(body)
  sel$promoter <- pro$score
  sel$transcriptBody <- body$score
  smallNumber <- max(c(1e-6, min(pro$score), min(body$score)))
  sel$log2meanCoverage <- log2(pro$score + smallNumber) + log2(body$score + smallNumber)
  sel$PT_score <- log2(pro$score + smallNumber) - log2(body$score + smallNumber)
  sel <- sel[order(sel$PT_score, decreasing = TRUE)]
  return(sel)
}

NFRscore <- function(obj, txs,
                     seqlev=intersect(seqlevels(obj), seqlevels(txs)),
                     nucleosomeSize=150, nucleosomeFreeSize=100){
  stopifnot(is(obj, "GAlignments"))
  stopifnot(is(txs, "GRanges"))
  obj <- as(obj, "GRanges")
  mcols(obj) <- NULL
  obj <- promoters(obj, upstream = 0, downstream = 1)
  cvg <- coverage(obj)
  cvg <- cvg[sapply(cvg, mean)>0]
  cvg <- cvg[names(cvg) %in% seqlev]
  seqlev <- seqlev[seqlev %in% names(cvg)]
  cvg <- cvg[seqlev]
  txs <- txs[seqnames(txs) %in% seqlev]
  txs <- unique(txs)
  sel <- promoters(txs, upstream = nucleosomeSize + floor(nucleosomeFreeSize/2),
                      downstream = nucleosomeSize + ceiling(nucleosomeFreeSize/2))

  n1.gr <- promoters(sel, upstream = 0, downstream = nucleosomeSize)
  n2.gr <- shift(n1.gr, shift = nucleosomeSize + nucleosomeFreeSize)
  n2.gr[strand(n2.gr)=="-"] <- shift(n1.gr[strand(n2.gr)=="-"],
                                     shift = -1 * (nucleosomeSize + nucleosomeFreeSize))
  nf.gr <- shift(n1.gr, shift = nucleosomeSize)
  width(nf.gr) <- nucleosomeFreeSize
  nf.gr[strand(nf.gr)=="-"] <- shift(n1.gr[strand(nf.gr)=="-"], shift = -1 * nucleosomeSize)
  start(nf.gr[strand(nf.gr)=="-"]) <- end(nf.gr[strand(nf.gr)=="-"]) - nucleosomeFreeSize + 1

  n1.gr$source <- "n1"
  n2.gr$source <- "n2"
  nf.gr$source <- "nf"
  n1.gr$oid <- seq_along(n1.gr)
  n2.gr$oid <- seq_along(n2.gr)
  nf.gr$oid <- seq_along(nf.gr)

  sel.gr <- c(n1.gr, nf.gr, n2.gr)
  sel.gr <- sel.gr[order(sel.gr$oid)]

  sel.gr <- split(sel.gr, seqnames(sel.gr))
  seqlev <- seqlev[seqlev %in% names(sel.gr)]
  sel.gr <- sel.gr[seqlev]
  cvg <- cvg[seqlev]
  vws <- Views(cvg, sel.gr)
  vms <- viewMeans(vws)
  sel.gr <- unlist(sel.gr)
  sel.gr$score <- unlist(vms)
  n1 <- sel.gr[sel.gr$source %in% "n1"]
  n2 <- sel.gr[sel.gr$source %in% "n2"]
  nf <- sel.gr[sel.gr$source %in% "nf"]
  stopifnot(identical(n1$oid, n2$oid))
  stopifnot(identical(n1$oid, nf$oid))
  sel <- sel[n1$oid]
  sel$n1 <- n1$score
  sel$nf <- nf$score
  sel$n2 <- n2$score
  smallNumber <- max(c(1e-6, min(nf$score), min(n1$score), min(n2$score)))
  sel$log2meanCoverage <- log2((3 * (n1$score + n2$score) + 2 * nf$score)/8 + smallNumber)
  sel$NFR_score <- log2(nf$score + smallNumber) + 1 - log2(n1$score + n2$score + smallNumber)
  sel <- sel[order(sel$NFR_score, decreasing = TRUE)]
  return(sel)
}


biobarplot <- function(x, xlab = "", ylab = "", col = NA) {
    sample.labels <- names(x)
    condition.labels <- colnames(x[[1]])

    if (is.na(col)) {
        col <- c("#666666", "#CCCCCC")
    }

    # prepare variables to save mean, sd, and p-values
    dfm <- dfs <- matrix(0, ncol = length(condition.labels), nrow = length(sample.labels))
    colnames(dfm) <- colnames(dfs) <- condition.labels
    rownames(dfm) <- rownames(dfs) <- sample.labels
    pvalues <- rep(NA, length(sample.labels))

    # calculate mean, sd, and p-values
    for (i in seq(x)) {
        dfm[i, ] <- apply(x[[i]], 2, mean, na.rm = TRUE)
        dfs[i, ] <- apply(x[[i]], 2, sd, na.rm = TRUE)
        x1 <- x[[i]][, 1]
        x2 <- x[[i]][, 2]
        pvalues[i] <- t.test(x1[!is.na(x1)], x2[!is.na(x2)])$p.value
    }

    # change data structure
    dfm <- t(dfm)
    dfs <- t(dfs)

    # calculate the y-coordinates for plotting *
    maxy <- max(dfm + dfs) * 1.1
    stepy <- max(dfm + dfs) * 0.1

    # bar chart
    bb <- barplot(dfm, beside = TRUE, ylim = c(0, maxy + 2 * stepy), col = col, border = col, xlab = xlab, ylab = ylab)

    # error bar
    arrows(bb, dfm - dfs, bb, dfm + dfs, code = 3, lwd = 1, angle = 90, length = 0.25 / length(sample.labels))

    # write *
    for (i in 1:length(pvalues)) {
        xi <- bb[, i]
        yi <- dfm[, i] + stepy * 1.5
        maxyi <- max(yi) + stepy
        if (pvalues[i] < 0.05) {
            lines(c(xi[1], xi[1], xi[2], xi[2]), c(yi[1], maxyi, maxyi, yi[2]))
            if (pvalues[i] < 0.01) {
                text((xi[1] + xi[2]) / 2, maxyi + stepy / 4, "**")
            } else if (pvalues[i] < 0.05) {
                text((xi[1] + xi[2]) / 2, maxyi + stepy / 4, "*")
            }
        }
    }

    # graph legend
    legend("topleft", legend = condition.labels, fill = col, col = col, border = col, box.lwd = 0, box.lty = 0)
}

PATHWAY_ENRICHMENT=function(project="NA",DEG_LIST=NA,GENE_LIST=NA){
 PATHWAY_LIST=read.table("~/archive_indiv/reference/IPA/180530_pathway_gene_list_2.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")

 FILTER_PATH=function(x){as.numeric_f(x[2])*as.numeric_f(x[3])!=0}
 DHYPER=function(x){1-phyper(as.numeric_f(x[2])-1, as.numeric_f(x[3]), as.numeric_f(x[4]), as.numeric_f(x[5]))}
 ODDS=function(x){as.numeric_f(x[2])*(as.numeric_f(x[3])+as.numeric_f(x[4])-as.numeric_f(x[5]))/(as.numeric_f(x[5])*(as.numeric_f(x[3])-as.numeric_f(x[2])))}

 dir.create_p("Enrichment_Analysis_res/full_res")
 #dir.create_p("Enrichment_Analysis_res/res")
 dir.create_p("Enrichment_Analysis_res/res_2")
 #dir.create_p("Enrichment_Analysis_res/res_3")

 pathway_gene_list=strsplit(PATHWAY_LIST$gene,",")

 RES=data.frame(PATHWAY=PATHWAY_LIST$pathway)
 RES$setdiff_GENE_PATH = RES$intersect_GENE_PATH = RES$intersect_DEG_PATH = NA

 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_DEG_PATH[kkk]=length(intersect(DEG_LIST,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_GENE_PATH[kkk]=length(intersect(GENE_LIST,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$setdiff_GENE_PATH[kkk]=length(setdiff(GENE_LIST,pathway_gene_list[[kkk]]))}
 RES$num_DEG=length(DEG_LIST)

 RES_2=RES[apply(RES,1,FILTER_PATH),]
 RES_2$Pvalue_hyper=apply(RES_2,1,DHYPER)
 RES_2$logPvalue=round(-log10(RES_2$Pvalue_hyper),digit=3)
 RES_2$Qvalue_BH=p.adjust(RES_2$Pvalue_hyper, method = "BH")
 RES_2$ODDS=apply(RES_2,1,ODDS)
 RES_2$RATIO=round(RES_2$intersect_DEG_PATH/RES_2$intersect_GENE_PATH,digits=3)
 RES_3=RES_2[order(RES_2$Pvalue_hyper),]

 write.table_FT(RES_3,paste0("Enrichment_Analysis_res/full_res/",today(),"_",project,"_DEG_Pathway_full.txt"))
 #write.table_FT(RES_3[RES_3$Qvalue_BH<0.05,c(1,6:10)],paste0("Enrichment_Analysis_res/res/",today(),"_",project,"_DEG_Pathway.txt"))
 write.table_FT(RES_3[(RES_3$Qvalue_BH<0.05)&(RES_3$intersect_GENE_PATH>5),c(1,6:10)],paste0("Enrichment_Analysis_res/res_2/",today(),"_",project,"_DEG_Pathway_filter2.txt"))
 #write.table_FT(RES_3[(RES_3$Pvalue_hyper<0.05)&(RES_3$intersect_GENE_PATH>1),c(1,6:10)],paste0("Enrichment_Analysis_res/res_3/",today(),"_",project,"_DEG_Pathway_filter3.txt"))
 }

PATHWAY_ENRICHMENT_2=function(project="NA",DEG_LIST=NA,GENE_LIST=NA){
 PATHWAY_LIST=read.table("~/archive_indiv/reference/IPA/180530_pathway_gene_list_2.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")

 FILTER_PATH=function(x){as.numeric_f(x[2])*as.numeric_f(x[3])!=0}
 DHYPER=function(x){1-phyper(as.numeric_f(x[2])-1, as.numeric_f(x[3]), as.numeric_f(x[4]), as.numeric_f(x[6]))}
 ODDS=function(x){as.numeric_f(x[2])*(as.numeric_f(x[3])+as.numeric_f(x[4])-as.numeric_f(x[6]))/(as.numeric_f(x[6])*(as.numeric_f(x[3])-as.numeric_f(x[2])))}

 pathway_gene_list=strsplit(PATHWAY_LIST$gene,",")

 RES=data.frame(PATHWAY=PATHWAY_LIST$pathway)
 RES$setdiff_GENE_PATH = RES$intersect_GENE_PATH = RES$intersect_DEG_PATH  = NA
 RES$intersect_DEG_name = NA

 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_DEG_PATH[kkk]=length(intersect(DEG_LIST,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_GENE_PATH[kkk]=length(intersect(GENE_LIST,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$setdiff_GENE_PATH[kkk]=length(setdiff(GENE_LIST,pathway_gene_list[[kkk]]))}
 for(kkk in 1:nrow(PATHWAY_LIST)){RES$intersect_DEG_name[kkk]=intersect(DEG_LIST,pathway_gene_list[[kkk]]) %>% paste0(.,collapse=",")}
 RES$num_DEG=length(DEG_LIST)

 RES_2=RES[apply(RES,1,FILTER_PATH),]
 RES_2$Pvalue_hyper=apply(RES_2,1,DHYPER)
 RES_2$logPvalue=round(-log10(RES_2$Pvalue_hyper),digit=3)
 RES_2$Qvalue_BH=p.adjust(RES_2$Pvalue_hyper, method = "BH")
 RES_2$ODDS=apply(RES_2,1,ODDS)
 RES_2$RATIO=round(RES_2$intersect_DEG_PATH/RES_2$intersect_GENE_PATH,digits=3)
 RES_3=RES_2[order(RES_2$Pvalue_hyper),]
 RES_3 = RES_3[,c("PATHWAY","intersect_DEG_PATH","intersect_GENE_PATH",
                  "setdiff_GENE_PATH","num_DEG",
                  "Pvalue_hyper","logPvalue","Qvalue_BH",
                  "ODDS","RATIO","intersect_DEG_name")]

 write.table_FT_2(RES_3,paste0("Enrichment_Analysis_res/full_res/",today(),"_",project,"_DEG_Pathway_full.txt"))
 write.table_FT_2(RES_3[(RES_3$Qvalue_BH<0.05)&(RES_3$intersect_GENE_PATH>5),c(1,6:10)],paste0("Enrichment_Analysis_res/res_2/",today(),"_",project,"_DEG_Pathway_filter2.txt"))
 }


is.blank <- function(x) {is.na(x) | x == ""}

distCor <- function(x) {as.dist(1-cor(x))}
zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
  if (scale=="row"){z <- t(scale(t(x)))}
  if (scale=="col"){z <- scale(x)}
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  hcl_row <- hclust(distCor(t(z)), method=method)
  hcl_col <- hclust(distCor(z), method=method)
  return(list(data=z, Rowv=as.dendrogram(hcl_row), Colv=as.dendrogram(hcl_col)))
  }

ggColorHue = function(n, l=65) {
              hues = seq(15, 375, length=n+1)
              hcl(h=hues, l=l, c=100)[1:n]}

file_list  = function(PATH_LIST){
              list      = data.frame(PATH=PATH_LIST)
              list$FILE = basename(list$PATH)
              print(list)
              }

bio_inst   = function(PACKAGE){
              source("https://bioconductor.org/biocLite.R")
              biocLite(PACKAGE)}
bio_inst_2   = function(PACKAGE){
                BiocManager::install(PACKAGE)}


make_list  = function(DIR,FILE,tag=NA){
              list = data.frame(PATH=list.files_n(DIR,FILE))
              list$FILE =basename(list$PATH)
              if(!is.na(tag)){colnames(list)=paste0(colnames(list),"_",tag)}
              list}

headd     = function(data.frame,num=5){
              if(ncol(data.frame)<num){col_num=ncol(data.frame)}else{col_num=num}
              if(nrow(data.frame)<num){row_num=nrow(data.frame)}else{row_num=num}
              data.frame[1:row_num,1:col_num]
             }

wilcox = function(d1,d2,NAME="no_name"){
                 suppressPackageStartupMessages(library(coin))
                 suppressPackageStartupMessages(library(DBI))
                 if(length(d1)*length(d2)!=0){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    tmp_res     =  wilcox_test(c(d1,d2) ~ factor(c(rep("d1",length(d1)),rep("d2",length(d2)))),distribution="exact")
                    up_down     = if(mean(d1)<mean(d2)){x=1;x}else{if(mean(d1)>mean(d2)){x=-1;x}else{x=0;x}}
                    Pvalue      = pvalue(tmp_res) %>% as.numeric_f()
                    Pvalue_SN   = SN_convert(Pvalue)
                    Pvalue_mark = p_convert(Pvalue)
                    tmp=data.frame(up_down=up_down,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(up_down=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                  }
                  rownames(tmp)=NAME
                  tmp
               }

wilcox_2 = function(d1,d2,NAME="no_name"){
                 suppressPackageStartupMessages(library(coin))
                 suppressPackageStartupMessages(library(DBI))
                 if((length(d1)*length(d2)!=0)&((unique(d1)%>%length()-1)+(unique(d2)%>%length()-1)!=0)){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    tmp_res     =  wilcox_test(c(d1,d2) ~ factor(c(rep("d1",length(d1)),rep("d2",length(d2)))),distribution="exact")
                    up_down     = if(mean(d1)<mean(d2)){x=1;x}else{if(mean(d1)>mean(d2)){x=-1;x}else{x=0;x}}
                    Pvalue      = pvalue(tmp_res) %>% as.numeric_f()
                    Pvalue_SN   = SN_convert(Pvalue,2,"E")
                    Pvalue_mark = p_convert(Pvalue)
                    tmp=data.frame(up_down=up_down,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(up_down=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                  }
                  rownames(tmp)=NAME
                  tmp
               }

student_t = function(d1,d2,NAME="no_name"){
                 if((length(d1)*length(d2)!=0)&((unique(d1)%>%length()-1)+(unique(d2)%>%length()-1)!=0)){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    up_down     = if(mean(d1)<mean(d2)){x=1;x}else{if(mean(d1)>mean(d2)){x=-1;x}else{x=0;x}}
                    Pvalue      = t.test(d1,d2, var.equal=F)$p.value %>% as.numeric_f() 
                    Pvalue_SN   = SN_convert(Pvalue,2,"E")
                    Pvalue_mark = p_convert(Pvalue)
                    tmp=data.frame(up_down=up_down,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(up_down=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                  }
                  rownames(tmp)=NAME
                  tmp
               }

##  wilcox = function(d1,d2,NAME="no_name"){
##                   suppressPackageStartupMessages(library(coin))
##                   if(length(d1)*length(d2)!=0){
##                      d1 = as.numeric(d1)
##                      d2 = as.numeric(d2)
##                      tmp_res     =  wilcox_test(c(d1,d2) ~ factor(c(rep("d1",length(d1)),rep("d2",length(d2)))),distribution="exact")
##                      up_down     = if(mean(d1)<mean(d2)){x=1;x}else{if(mean(d1)>mean(d2)){x=-1;x}else{x=0;x}}
##                      Pvalue      = pvalue(tmp_res)
##                      Pvalue_SN   = SN_convert(Pvalue)
##                      Pvalue_mark = p_convert(Pvalue)
##                      tmp=data.frame(up_down=up_down,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
##                      tmp=data.frame(up_down=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
##                    }
##                    rownames(tmp)=NAME
##                    tmp
##                 }

## wilcox_2 = function(d1,d2,NAME="no_name"){
##                  suppressPackageStartupMessages(library(coin))
##                  if(length(d1)*length(d2)!=0){
##                     d1 = as.numeric(d1)
##                     d2 = as.numeric(d2)
##                     tmp_res     =  wilcox_test(c(d1,d2) ~ factor(c(rep("d1",length(d1)),rep("d2",length(d2)))),distribution="exact")
##                     up_down     = if(mean(d1)<mean(d2)){x=1;x}else{if(mean(d1)>mean(d2)){x=-1;x}else{x=0;x}}
##                     Pvalue      = pvalue(tmp_res) %>% as.numeric_f()
##                     Pvalue_SN   = SN_convert(Pvalue,2,"e")
##                     Pvalue_mark = p_convert(Pvalue)
##                     tmp=data.frame(up_down=up_down,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
##                     tmp=data.frame(up_down=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
# #                  }
##                   rownames(tmp)=NAME
##                   tmp
##                }

pearson = function(d1,d2,NAME="no_name"){
                   if(length(d1)*length(d2)!=0){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    Pearson     = cor.test(d1,d2,method="pearson")
                    Pvalue      = Pearson$p.value
                    cor         = Pearson$estimate %>% SN_convert(.,2,"xE")
                    Pvalue_mark = p_convert(Pvalue)
                    Pvalue_SN   = SN_convert(Pvalue,2,"xE")
                    tmp=data.frame(cor=cor,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(cor=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                   }
                  rownames(tmp)=NAME
                  tmp
                 }

pearson_2 = function(d1,d2,NAME="no_name"){
                   if(length(d1)*length(d2)!=0){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    Pearson     = cor.test(d1,d2,method="pearson")
                    Pvalue      = Pearson$p.value
                    cor         = Pearson$estimate %>% SN_convert(.,2,"e")
                    Pvalue_mark = p_convert(Pvalue)
                    Pvalue_SN   = SN_convert(Pvalue,2,"e")
                    tmp=data.frame(cor=cor,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(cor=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                   }
                  rownames(tmp)=NAME
                  tmp
                 }

spearman = function(d1,d2,NAME="no_name"){
                   if(length(d1)*length(d2)!=0){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    Pearson     = cor.test(d1,d2,method="spearman")
                    Pvalue      = Pearson$p.value
                    rho         = Pearson$estimate %>% SN_convert(.,2,"xE")
                    Pvalue_mark = p_convert(Pvalue)
                    Pvalue_SN   = SN_convert(Pvalue,2,"xE")
                    tmp=data.frame(rho=rho,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(rho=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                   }
                   rownames(tmp)=NAME
                   tmp
                 }

spearman_2 = function(d1,d2,NAME="no_name"){
                   if(length(d1)*length(d2)!=0){
                    d1 = as.numeric(d1)
                    d2 = as.numeric(d2)
                    Pearson     = suppressWarnings(cor.test(d1,d2,method="spearman"))
                    Pvalue      = Pearson$p.value
                    rho         = Pearson$estimate %>% SN_convert(.,2,"e")
                    Pvalue_mark = p_convert(Pvalue)
                    Pvalue_SN   = SN_convert(Pvalue,2,"e")
                    tmp=data.frame(rho=rho,Pvalue=Pvalue,Pvalue_SN=Pvalue_SN,Pvalue_mark=Pvalue_mark)}else{
                    tmp=data.frame(rho=NA,Pvalue=1,Pvalue_SN=1,Pvalue_mark="")
                   }
                   rownames(tmp)=NAME
                   tmp
                 }


geom_tile_plot =function(data,colname="COLNAME",rowname="ROWNAME",Title=NA,
                         palette="Pastel1",color=c("#ff2800","#faf500","#35a16b","#0041ff","#66ccff","#ff99a0","#ff9900","#9a0079","#c7b2de","#663300","#7f878f","#000000","#ffffff"),
                         out_dir=NA,color_mode="palette",
                         color_bar=NA,color_bar_name="COLOR_BAR",color_bar_levels=NA,color_bar_pallete="Pastel1",color_bar_mode="palette",
                         color_bar_color=c("#ff2800","#faf500","#35a16b","#0041ff","#66ccff","#ff99a0","#ff9900","#9a0079","#c7b2de","#663300","#7f878f","#000000","#ffffff"),
                         angle_col=60,text_col=10,text_row=1,pdf_h=5,pdf_w=10,p1_h=0.8,p2_h=0.1){
                           suppressPackageStartupMessages(library(tidyr))
                           suppressPackageStartupMessages(library(ggplot2))

                           colnames(data) = gsub("-",".",colnames(data))
                           rownames(data) = gsub("-",".",rownames(data))

                           df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)
                           p= ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                               geom_tile()+
                               theme_bw() +
                               scale_x_discrete(limits=unique(df$colnames)) +
                               scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
                               theme(axis.text.x=element_text(colour="black",angle = angle_col,hjust=1,vjust=1,size=text_col,face="bold", family = "Helvetica"),
                                     axis.text.y=element_text(colour="black",size=text_row,face="bold", family = "Helvetica"),
                                     axis.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
                                     strip.text=element_text(size=8,face="bold", family = "Helvetica"),
                                     axis.line = element_line(colour = "black"),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                               labs(x=colname,y=rowname)

                           if(color_mode=="palette"){p = p+scale_fill_brewer(palette =palette)}else{
                                                     p = p+scale_fill_manual(values = color)  }
                           if(length(color_bar)==1){p=p+labs(title=paste0(today()," ",Title))}

                           if(length(color_bar)!=1){
                            df2 = data.frame(colnames=unique(df$colnames),BAR=color_bar)
                            if(is.na(color_bar_levels)){df2$BAR = factor(df2$BAR)}else{df2$BAR = factor(df2$BAR,levels=color_bar_levels)}

                            p2= ggplot(df2, aes(x=colnames,y=1,fill=color_bar))+
                                 geom_tile()+
                                 scale_y_continuous(expand=c(0,0)) +
                                 scale_x_discrete(limits=df2$colnames) +
                                 theme(axis.title.x=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),
                                       axis.text.x=element_blank(),axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
                                       plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
                                 labs(title=paste0(today()," ",Title),y=color_bar_name)
                            if(color_bar_mode=="palette"){p2=p2+scale_fill_brewer(color_bar_pallete =palette)}else{p2=p2+scale_fill_manual(values = color_bar_color)}
                           }

                           if(!is.na(out_dir)){dir.create_p(out_dir)}
                           if(is.na(out_dir)){out_dir=getwd()}

                           if(length(color_bar)==1){
                            pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                             plot(p)
                            dev.off()
                           }

                           if(length(color_bar)!=1){
                              suppressPackageStartupMessages(library(gridExtra))
                              gp1 = ggplotGrob(p)
                              gp2 = ggplotGrob(p2)

                              maxWidth = grid::unit.pmax(gp1$widths, gp2$widths)
                              gp1$widths <- gp2$widths <- as.list(maxWidth)

                              pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                               grid.arrange(gp2,gp1, ncol=1,heights=c(p2_h,p1_h))
                              dev.off()
                           }

                         }

geom_tile_plot_2_old =function(data,colname="COLNAME",rowname="ROWNAME",Title=NA,out_dir=NA,
                           angle_col=60,text_col=10,text_row=1,pdf_h=5,pdf_w=10){
                           df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)
                           p= ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                               geom_tile()+
                               scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0)+
                               theme_bw() +
                               scale_x_discrete(limits=colnames(data)) +
                               scale_y_discrete(limits=rownames(data)%>%rev()) +
                               theme(axis.text.x=element_text(colour="black",angle = angle_col,hjust=1,vjust=1,size=text_col,face="bold", family = "Helvetica"),
                                     axis.text.y=element_text(colour="black",size=text_row,face="bold", family = "Helvetica"),
                                     axis.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
                                     strip.text=element_text(size=8),
                                     axis.line = element_line(colour = "black"),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                               labs(title=paste0(today()," ",Title),x=colname,y=rowname)

                           if(is.na(out_dir)){out_dir=getwd()}
                           pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                            plot(p)
                           dev.off()
                         }

geom_tile_plot_2 =function(data,colname="COLNAME",rowname="ROWNAME",Title=NA,out_dir=NA,
                            color_bar=NA,color_bar_name="COLOR_BAR",color_bar_levels=NA,color_bar_pallete="Pastel1",color_bar_mode="palette",
                            color_bar_color=c("#ff2800","#faf500","#35a16b","#0041ff","#66ccff","#ff99a0","#ff9900","#9a0079","#c7b2de","#663300","#7f878f","#000000","#ffffff"),
                            angle_col=60,text_col=10,text_row=1,pdf_h=5,pdf_w=10,p1_h=0.8,p2_h=0.1){
                           suppressPackageStartupMessages(library(tidyr))
                           suppressPackageStartupMessages(library(ggplot2))

                           colnames(data) = gsub("-",".",colnames(data))
                           rownames(data) = gsub("-",".",rownames(data))

                           df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)
                           p= ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                               geom_tile()+
                               scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0)+
                               theme_bw() +
                               scale_x_discrete(limits=unique(df$colnames)) +
                               scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
                               theme(axis.text.x=element_text(colour="black",angle = angle_col,hjust=1,vjust=1,size=text_col,face="bold", family = "Helvetica"),
                                     axis.text.y=element_text(colour="black",size=text_row,face="bold", family = "Helvetica"),
                                     axis.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
                                     strip.text=element_text(size=8),
                                     axis.line = element_line(colour = "black"),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
                               labs(x=colname,y=rowname)

                           if(length(color_bar)==1){p=p+labs(title=Title)}

                           if(length(color_bar)!=1){
                            df2 = data.frame(colnames=unique(df$colnames),BAR=color_bar)
                            if(length(color_bar_levels)==1){df2$BAR = factor(df2$BAR)}else{df2$BAR = factor(df2$BAR,levels=color_bar_levels)}

                            if(color_bar_mode=="palette"){
                             p2= ggplot(df2, aes(x=colnames,y=1,fill=BAR))+
                                 geom_tile()+scale_fill_brewer(palette=color_bar_pallete)+
                                 scale_y_continuous(expand=c(0,0)) +
                                 scale_x_discrete(limits=df2$colnames) +
                                 theme(axis.title.x=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),
                                       axis.text.x=element_blank(),axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
                                       plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
                                 labs(title=paste0(today()," ",Title),y=color_bar_name)}else{
                                  p2= ggplot(df2, aes(x=colnames,y=1,fill=BAR))+
                                 geom_tile()+scale_fill_manual(values = color_bar_color)+
                                 scale_y_continuous(expand=c(0,0)) +
                                 scale_x_discrete(limits=df2$colnames) +
                                 theme(axis.title.x=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),
                                       axis.text.x=element_blank(),axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
                                       plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
                                 labs(title=paste0(today()," ",Title),y=color_bar_name)

                                 }
                               }

                           if(!is.na(out_dir)){dir.create_p(out_dir)}
                           if(is.na(out_dir)){out_dir=getwd()}

                           if(length(color_bar)==1){
                            pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                             plot(p)
                            dev.off()
                           }

                           if(length(color_bar)!=1){
                              suppressPackageStartupMessages(library(gridExtra))
                              gp1 = ggplotGrob(p)
                              gp2 = ggplotGrob(p2)

                              maxWidth = grid::unit.pmax(gp1$widths, gp2$widths)
                              gp1$widths <- gp2$widths <- as.list(maxWidth)

                              pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                               grid.arrange(gp2,gp1, ncol=1,heights=c(p2_h,p1_h))
                              dev.off()
                           }
                         }
## colnameの先頭文字に数字は禁則


geom_tile_plot_3 =function(data,colname="COLNAME",rowname="ROWNAME",Title=NA,out_dir=NA,
                           angle_col=60,text_col=10,text_row=1,pdf_h=5,pdf_w=10){

                           colnames(data) = gsub("-",".",colnames(data))
                           rownames(data) = gsub("-",".",rownames(data))

                           df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)
                           p= ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                               geom_tile()+
                               scale_fill_gradient(low = "white", high = "red")+
                               theme_bw() +
                               scale_x_discrete(limits=unique(df$colnames)) +
                               scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
                               theme(axis.text.x=element_text(colour="black",angle = angle_col,hjust=1,vjust=1,size=text_col,face="bold", family = "Helvetica"),
                                     axis.text.y=element_text(colour="black",size=text_row,face="bold", family = "Helvetica"),
                                     axis.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
                                     strip.text=element_text(size=8,face="bold", family = "Helvetica"),
                                     axis.line = element_line(colour = "black"),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                               labs(title=paste0(today()," ",Title),x=colname,y=rowname)

                           if(!is.na(out_dir)){dir.create_p(out_dir)}
                           if(is.na(out_dir)){out_dir=getwd()}

                           pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                            plot(p)
                           dev.off()
                         }

geom_tile_plot_FT =function(data,colname="COLNAME",rowname="ROWNAME",Title=NA,palette="Pastel1",out_dir=NA,
                            angle_col=60,text_col=10,text_row=10,pdf_h=5,pdf_w=10){
                           suppressPackageStartupMessages(library(tidyr))
                           suppressPackageStartupMessages(library(ggplot2))

                           colnames(data) = gsub("-",".",colnames(data))
                           rownames(data) = gsub("-",".",rownames(data))

                           df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)
                           p= ggplot(df,aes(x=colnames,y=rownames,fill=score))+
                               geom_tile()+
                               scale_fill_manual(values =c("white","red"))+
                               theme_bw() +
                               scale_x_discrete(limits=unique(df$colnames)) +
                               scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
                               theme(axis.text.x=element_text(colour="black",angle = angle_col,hjust=1,vjust=1,size=text_col,face="bold", family = "Helvetica"),
                                     axis.text.y=element_text(colour="black",size=text_row,face="bold", family = "Helvetica"),
                                     axis.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
                                     strip.text=element_text(size=8,face="bold", family = "Helvetica"),
                                     axis.line = element_line(colour = "black"),
                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none") +
                               labs(title=paste0(today()," ",Title),x=colname,y=rowname)

                           if(!is.na(out_dir)){dir.create_p(out_dir)}
                           if(is.na(out_dir)){out_dir=getwd()}

                           pdf(paste0(out_dir,"/",today(),"_",Title,"_geom_tile.pdf"),h=pdf_h,w=pdf_w)
                            plot(p)
                           dev.off()
                         }

hclust_col_order=function(DATA,METHOD="complete"){
  suppressPackageStartupMessages(library(ggdendro))
  rd     = t(DATA) %>% dist()
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  col_order = as.character(ddata$labels$label)
  col_order
}

hclust_col_order_ward=function(DATA,METHOD="ward.D2"){
  suppressPackageStartupMessages(library(ggdendro))
  rd     = t(DATA) %>% dist()
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  col_order = as.character(ddata$labels$label)
  col_order
}

hclust_row_order=function(DATA,METHOD="complete"){
  suppressPackageStartupMessages(library(ggdendro))
  rd     = dist(DATA)
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  row_order = as.character(ddata$labels$label)
  row_order
}

hclust_row_order_ward=function(DATA,METHOD="ward.D2"){
  suppressPackageStartupMessages(library(ggdendro))
  rd     = dist(DATA)
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  row_order = as.character(ddata$labels$label)
  row_order
}

exp_filter=function(DATA,exp=10,p=0.1){apply(DATA,1,function(x){length(x[x>=exp])>=length(x)*p})}

Module_reconstruction = function(res,min_match=0,min_module=5){
  for(iii in 1:ncol(res)){
    res[is.na(res)] = "ZZZ_NoData"
    res[,iii] = factor(res[,iii]) %>% as.numeric()
  }

  res_full = matrix(NA,nrow=nrow(res),ncol=nrow(res)) %>% as.data.frame()
  colnames(res_full) = rownames(res_full) = rownames(res)

  for(jjj in 1:nrow(res)){
    hoge     = rep(res[jjj,],nrow(res)) %>% as.numeric_f() %>% matrix(nrow(res),ncol(res),byrow=T)
    hogera   = res-hoge
    hogera2  = hogera==0
    res_full[jjj,] = rowSums(hogera2) %>% as.numeric_f()
    print(paste0(jjj," / ",nrow(res)))
  }

  #********************************************
  colnames(res_full) = c(1:ncol(res_full))
  rownames(res_full) = c(1:nrow(res_full))

  res_full = data.frame(Gene1=rownames(res_full),res_full)
  colnames(res_full) = gsub("^X","",colnames(res_full))

  res_full = gather(res_full,key=Gene2,value=match_num,-Gene1)
  res_full$Gene1 = as.numeric(res_full$Gene1)
  res_full$Gene2 = as.numeric(res_full$Gene2)

  res_full = res_full[res_full$Gene1<res_full$Gene2,]

  #********************************************

  #---------------------------------------------------------#
  ## 取り急ぎ、moduleは0個以上のsubsetで合致したものに限定してみる (==限定しない)
  #---------------------------------------------------------#

  res_full_part = res_full[res_full$match_num>=min_match,]
  dim(res_full_part)
  unique(c(res_full_part$Gene1,res_full_part$Gene2)) %>% sort() %>% length() ## 139

  #-----------------------------------
  for(iii in c(max(res_full_part$match_num):0)){
    res_tmp = res_full_part[res_full_part$match_num>=iii,]
    res_tmp = res_tmp[order(res_tmp$Gene1,res_tmp$Gene2),]
    if(iii!=max(res_full_part$match_num)){res_tmp = res_tmp[(!is.element(res_tmp$Gene1,Module_omit_gene))&(!is.element(res_tmp$Gene2,Module_omit_gene)),]}
    if(nrow(res_tmp)==0){print(paste0(iii," was 0 gene"))}
    if(nrow(res_tmp)==0){break}

    tmp_gene_list_1 = res_tmp$Gene1 %>% unique() %>% as.numeric_f() %>% sort()
    tmp_gene_list_2 = c(res_tmp$Gene1,res_tmp$Gene2) %>% unique() %>% as.numeric_f() %>% sort()

    Module_list = list()

    for(kkk in 1:length(tmp_gene_list_1)){
      target_tmp         = tmp_gene_list_1[kkk]
      tmp_row_num        = (target_tmp==res_tmp$Gene1)|(target_tmp==res_tmp$Gene2)
      gene_module_tmp    = res_tmp[tmp_row_num,c("Gene1","Gene2")] %>% unlist() %>% as.numeric_f() %>% unique() %>% sort()
      if(min(gene_module_tmp)==target_tmp){
        Module_list[[kkk]]=gene_module_tmp}else{
          lll=sapply(Module_list,function(x){is.element(min(gene_module_tmp),x)}) %>% grep("TRUE",.) %>% min()
          Module_list[[lll]]=union(Module_list[[lll]],gene_module_tmp) %>% sort() %>% unique()
        }
      }
    list_tmp      = data.frame(module_num=c(1:length(Module_list)),module_size=sapply(Module_list,length))
    Module_list_2 = Module_list[list_tmp$module_num[list_tmp$module_size>0]]

    Module_final =list()

    if(length(unlist(Module_list_2))==length(tmp_gene_list_2)){Module_final=Module_list_2}else{
     tmp_res           = matrix(NA,nrow=length(Module_list_2),ncol=length(tmp_gene_list_2))
     colnames(tmp_res) = paste0("gene_",tmp_gene_list_2)
     rownames(tmp_res) = paste0("Module_",c(1:nrow(tmp_res)))
     for(nnn in 1:nrow(tmp_res)){
       tmp_res[nnn,] = is.element(tmp_gene_list_2,Module_list_2[[nnn]])
     }

    for(mmm in 1:ncol(tmp_res)){
      if(grep("TRUE",tmp_res[,mmm]) %>% length() >1){
       module_num_tmp              = grep("TRUE",tmp_res[,mmm])
       hoge                        = tmp_res[module_num_tmp,]
       hogera                      = apply(hoge,2,function(x){any(x=="TRUE")})
       tmp_res[module_num_tmp[1],] = hogera
       tmp_res[module_num_tmp[-1],] = FALSE
     }
    }

    tmp_res = tmp_res[apply(tmp_res,1,function(x){any(x=="TRUE")}),,drop=F]

    for(ppp in 1:nrow(tmp_res)){
      Module_final[[ppp]] = tmp_gene_list_2[tmp_res[ppp,]]}
    }

    list_tmp         = data.frame(module_num=c(1:length(Module_final)),module_size=sapply(Module_final,length))
    if(iii==max(res_full_part$match_num)){Module_final_2   = Module_final[list_tmp$module_num[list_tmp$module_size>=5]]}else{
     Module_final_2   = c(Module_final_2,Module_final[list_tmp$module_num[list_tmp$module_size>=min_module]])
    }

    Module_omit_gene = unlist(Module_final_2) %>% sort() %>% unique()

    print(paste0(iii," subset_Module_omit_gene ",length(Module_omit_gene)," genes"))
    print(paste0(iii," subset_Module_omit_gene ",length(Module_final_2),"modules"))
  }

  gene_corresp = data.frame(Gene=rownames(res))
  gene_corresp$module_num = NA
  for(iii in 1:length(Module_final_2)){
    gene_corresp$module_num[Module_final_2[[iii]]] = paste0("Module_",formatC(iii,width=2,flag="0"))
  }

  gene_corresp
}

grid_extra_2 = function(p1=p1,p2=p2,
                        p1_h=.2,p2_h=.8,
                        h_pdf=5,w_pdf=10,out_file="grid_extra.pdf"){
                          suppressPackageStartupMessages(library(gridExtra))
                          gp1 = ggplotGrob(p1)
                          gp2 = ggplotGrob(p2)

                          maxWidth = grid::unit.pmax(gp1$widths, gp2$widths)
                          gp1$widths <- gp2$widths <- as.list(maxWidth)

                          dir.create_p(dirname(out_file))

                          pdf(out_file,h=h_pdf,w=w_pdf)
                           grid.arrange(gp1,gp2, ncol=1,heights=c(p1_h,p2_h))
                          dev.off()
                        }

grid_extra_3 = function(p1=p1,p2=p2,p3=p3,
                        p1_h=.2,p2_h=.1,p3_h=.8,
                        h_pdf=5,w_pdf=10,out_file="grid_extra.pdf"){
                          suppressPackageStartupMessages(library(gridExtra))
                          gp1 = ggplotGrob(p1)
                          gp2 = ggplotGrob(p2)
                          gp3 = ggplotGrob(p3)

                          maxWidth = grid::unit.pmax(gp1$widths, gp2$widths, gp3$widths)
                          gp1$widths <- gp2$widths <- gp3$widths <- as.list(maxWidth)

                          dir.create_p(dirname(out_file))

                          pdf(out_file,h=h_pdf,w=w_pdf)
                           grid.arrange(gp1,gp2,gp3, ncol=1,heights=c(p1_h,p2_h,p3_h))
                          dev.off()
                        }

grid_extra_4 = function(p1=p1,p2=p2,p3=p3,p4=p4,
                        p1_h=.2,p2_h=.1,p3_h=.1,p4_h=.8,
                        h_pdf=5,w_pdf=10,out_file="grid_extra.pdf"){
                          suppressPackageStartupMessages(library(gridExtra))
                          gp1 = ggplotGrob(p1)
                          gp2 = ggplotGrob(p2)
                          gp3 = ggplotGrob(p3)
                          gp4 = ggplotGrob(p4)

                          maxWidth = grid::unit.pmax(gp1$widths, gp2$widths, gp3$widths, gp4$widths)
                          gp1$widths <- gp2$widths <- gp3$widths <- gp4$widths <- as.list(maxWidth)

                          dir.create_p(dirname(out_file))

                          pdf(out_file,h=h_pdf,w=w_pdf)
                           grid.arrange(gp1,gp2,gp3,gp4, ncol=1,heights=c(p1_h,p2_h,p3_h,p4_h))
                          dev.off()
                        }

grid_extra_5 = function(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,
                        p1_h=.2,p2_h=.1,p3_h=.1,p4_h=.1,p5_h=.8,
                        h_pdf=5,w_pdf=10,out_file="grid_extra.pdf"){
                          suppressPackageStartupMessages(library(gridExtra))
                          gp1 = ggplotGrob(p1)
                          gp2 = ggplotGrob(p2)
                          gp3 = ggplotGrob(p3)
                          gp4 = ggplotGrob(p4)
                          gp5 = ggplotGrob(p5)

                          maxWidth = grid::unit.pmax(gp1$widths, gp2$widths, gp3$widths, gp4$widths, gp5$widths)
                          gp1$widths <- gp2$widths <- gp3$widths <- gp4$widths <- gp5$widths <- as.list(maxWidth)

                          dir.create_p(dirname(out_file))

                          pdf(out_file,h=h_pdf,w=w_pdf)
                           grid.arrange(gp1,gp2,gp3,gp4,gp5, ncol=1,heights=c(p1_h,p2_h,p3_h,p4_h,p5_h))
                          dev.off()
                        }

grid_extra_6 = function(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,
                        p1_h=.2,p2_h=.1,p3_h=.1,p4_h=.1,p5_h=.1,p6_h=.8,
                        h_pdf=5,w_pdf=10,out_file="grid_extra.pdf"){
                          suppressPackageStartupMessages(library(gridExtra))
                          gp1 = ggplotGrob(p1)
                          gp2 = ggplotGrob(p2)
                          gp3 = ggplotGrob(p3)
                          gp4 = ggplotGrob(p4)
                          gp5 = ggplotGrob(p5)
                          gp6 = ggplotGrob(p6)

                          maxWidth = grid::unit.pmax(gp1$widths, gp2$widths, gp3$widths, gp4$widths, gp5$widths, gp6$widths)
                          gp1$widths <- gp2$widths <- gp3$widths <- gp4$widths <- gp5$widths <- gp6$widths <- as.list(maxWidth)

                          dir.create_p(dirname(out_file))

                          pdf(out_file,h=h_pdf,w=w_pdf)
                           grid.arrange(gp1,gp2,gp3,gp4,gp5,gp6, ncol=1,heights=c(p1_h,p2_h,p3_h,p4_h,p5_h,p6_h))
                          dev.off()
                        }


find_order = function(target,list){is.element(list,target)%>%grep("TRUE",.)}

sign_pick_tmp  = function(x){if(is.na(x)){x=0}else{if(x>0){x=1}else{if(x<0){x=-1}}};x}
sign_pick      = function(x){sapply(x,sign_pick_tmp)%>% as.numeric_f()}
sign_pick_2    = function(data.frame){res = matrix(NA,ncol=ncol(data.frame),nrow=nrow(data.frame)) %>% as.data.frame()
                                      colnames(res) = colnames(data.frame)
                                      rownames(res) = rownames(data.frame)

                                      for(iii in 1:nrow(data.frame)){
                                        res[iii,] = sign_pick(data.frame[iii,])
                                        }
                                      res
                                      }

circos_plot = function(data_1,data_2,NAME="tmp"){
                       suppressPackageStartupMessages(library(circlize))
                       suppressPackageStartupMessages(library(RColorBrewer))

                       df   = table(data_1,data_2)
                       df_2 = as.data.frame(df) %>% spread(key=data_2,value=Freq)
                       df_3 = df_2[,-1]
                       rownames(df_3) = df_2[,1]
                       colnames(df_3) = colnames(df_3)

                       col_list = colorRampPalette(c("red", "yellow", "green", "blue"))

                       pdf(paste0(today(),"_",NAME,"_circos_plot.pdf"))
                        chordDiagram(as.matrix(df_3),
                                     transparency = 0.5, #透明度の調整
                                     column.col   = col_list(ncol(df)), #行の色を設定、データ行で色数を指定
                                     row.col      = col_list(nrow(df)), #列の色を設定、データ列で色数を指定
                                     grid.col     = c(col_list(nrow(df)), col_list(ncol(df))) #ラベルの色を設定
                                     )
                        title(NAME)
                       dev.off()
                   }

# Fisherの正確確率検定
FISHER_p = function(class1,class2){
            tmp    = table(class1,class2) %>% as.data.frame()
            x      = matrix(tmp$Freq, ncol=2, byrow=T)
            FISHER = fisher.test(x)$p.value
            FISHER
            }


## 低発現変動遺伝子のfilter→TMM→logCPM
LowExp_TMM_logCPM = function(DATA){
                      suppressPackageStartupMessages(library(edgeR))

                      DATA_f  = DATA[exp_filter(DATA),]
                      DGE     = DGEList(counts=DATA_f, genes=rownames(DATA_f))
                      DGE_fn  = calcNormFactors( DGE ,method="TMM")
                      logCPM  = log2(cpm(DGE_fn)+1)
                      logCPM
                    }

## 相補的配列を得る
primer_res=function(FILE=FILE,OUT_FILE=OUT_FILE){
                   data=read.table_FT(FILE)
                   data$after = NA

                   for(iii in 1:nrow(data)){
                    tmp_list=data.frame(before=strsplit(data$ARRAY[iii],"")[[1]])
                    corresp =data.frame(before=c("A","G","C","T"),after=c("T","C","G","A"))

                    tmp_list_2 = left_join(tmp_list,corresp,by="before")
                    data$after[iii]=rev(tmp_list_2$after) %>% paste0(.,collapse="")
                   }

                   write.table_FT(data,OUT_FILE)}

######## clusterのelbow解析に使用
sqr_edist <- function(x, y) {sum((x-y)^2)}

wss.cluster <- function(clustermat) {
                 c0 <- apply(clustermat, 2, FUN=mean)
                 sum(apply(clustermat, 1, FUN=function(row){sqr_edist(row,c0)}))
                 }

wss.total <- function(dmatrix, labels) {
              wsstot <- 0
              k <- length(unique(labels))
              for(i in 1:k){
                  wsstot <- wsstot + wss.cluster(subset(dmatrix, labels==i))
               }
              wsstot
             }

elbow <- function(dmatrix, kmax) {
          npts <- dim(dmatrix)[1]
          wss <- numeric(kmax)
          wss[1] <- (npts-1)*sum(apply(dmatrix, 2, var))

          for(k in 2:kmax) {
              d <- dist(dmatrix, method="euclidean")
              pfit <- hclust(d, method="ward.D2")
              labels <- cutree(pfit, k=k)
              wss[k] <- wss.total(dmatrix, labels)
          }
          list(wss = wss)
         }
######## clusterのelbow解析に使用 ここまで #####

hum_to_mice_conv=function(Gene_set){paste0(substr(Gene_set,1,1),substr(Gene_set,2,nchar(Gene_set))%>%tolower())}


###########################################

## ref) https://to-kei.net/r-beginner/step-w/
## https://to-kei.net/r-beginner/r-3/
## http://cogpsy.educ.kyoto-u.ac.jp/personal/Kusumi/datasem13/shinya.pdf

step_mra_fun = function(dat,fn){
  colnames(dat)[1] = "y"
  ans   = lm(dat$y~., data=dat)
  ans   = step(ans)
  s.ans = summary(ans)
  coe   = s.ans$coefficient
  N     = nrow(dat)
  aic   = AIC(ans)
  conf  = confint(ans,levels=(0.95))
  R2    = summary(ans)$r.squared
  result= cbind(coe,conf,R2,aic,N)
  colnames(result)[5:6] = c("RR95%CI.low","RR95%CI.up")
  e2    = deviance(ans)
  MSE   = e2/N

  if(nrow(result)>=2){result[2:nrow(result),7:9] <- ""}

  write.table(cbind("E2",e2,"MSE",MSE),fn,append=T,quote=F,sep=",",row.names=F,col.names=F)
  write.table(matrix(c("",colnames(result)),nrow=1),fn,append=T,quote=F,sep=",",row.names=F,col.names=F)
  write.table(result,fn,append=T,quote=F,sep=",",row.names=T,col.names=F)
  write.table("",fn,append=T,quote=F,sep=",",row.names=F,col.names=F)
}

## convertMouseGeneList <- function(x){require("biomaRt")
##                                     human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
##                                     mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
## 
##                                     genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
##                                     humanx <- unique(genesV2[, 2])
## 
##                                     return(humanx)
##                                   }
convertMouseGeneList <- function(x){require("biomaRt")
                                    mart      = useMart(biomart = "ENSEMBL_MART_ENSEMBL")
                                    m_ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)
                                    h_ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

                                    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = m_ensembl, attributesL = c("hgnc_symbol"), martL = h_ensembl, uniqueRows=T)
                                    humanx <- unique(genesV2[, "HGNC.symbol"])

                                    return(humanx)
                                  }


# Basic function to convert human to mouse gene names
##    convertHumanGeneList <- function(x){require("biomaRt")
##                                        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
##                                        mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
##    
##                                        genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
##                                        humanx <- unique(genesV2[, 2])
##    
##                                        return(humanx)
##                                        }
convertHumanGeneList <- function(x){require("biomaRt")
                                    mart      = useMart(biomart = "ENSEMBL_MART_ENSEMBL")
                                    m_ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)
                                    h_ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

                                    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = h_ensembl, attributesL = c("mgi_symbol"), martL = m_ensembl, uniqueRows=T)
                                    mousex <- unique(genesV2[, "MGI.symbol"])

                                    return(mousex)
                                    }



###########################################
list_summarize = function(list){TARGET_full = unlist(list) %>% sort() %>% unique()
                                res = matrix(NA,ncol=length(list),nrow=length(TARGET_full)) %>% as.data.frame()
                                colnames(res) = names(list)
                                rownames(res) = TARGET_full
                                for(iii in 1:length(list)){res[,iii]=is.element(TARGET_full,list[[iii]])}
                                res
                               }

###########################################
table_freq       = function(data.frame,name="factor"){res_tmp = table(data.frame) %>% as.data.frame() %>% filter(Freq!=0);res_tmp$data.frame=as.character(res_tmp$data.frame);colnames(res_tmp)[1]=name;res_tmp}
as.data.frame_2  = function(table){tmp   = as.matrix(table)
                                   tmp_2 = matrix(NA,ncol=ncol(tmp),nrow=nrow(tmp)) %>% as.data.frame()
                                   colnames(tmp_2) = colnames(table)
                                   rownames(tmp_2) = rownames(table)
                                   for(iii in 1:ncol(tmp)){tmp_2[,iii]=tmp[,iii]}
                                   tmp_2
                                  }

table_to_dataframe = function(table){write.table(table,"tmptmptmptmptmptmp.txt")
                                     hoge=read.table("tmptmptmptmptmptmp.txt",header=T,row.names=1,stringsAsFactors=F)
                                     file.remove("tmptmptmptmptmptmp.txt")
                                     hoge}
table_freq_2 = function(df1,df2){table(df1,df2)%>%table_to_dataframe()}

###########################################
biobarplot <- function(x, xlab = "", ylab = "", col = NA) {
    sample.labels <- names(x)
    condition.labels <- colnames(x[[1]])

    if (is.na(col)) {
        col <- c("black", "white")
    }
    
    # prepare variables to save mean, sd, and p-values
    dfm <- dfs <- matrix(0, ncol = length(condition.labels), nrow = length(sample.labels))
    colnames(dfm) <- colnames(dfs) <- condition.labels
    rownames(dfm) <- rownames(dfs) <- sample.labels
    pvalues <- rep(NA, length(sample.labels))
    
    # calculate mean, sd, and p-values
    for (i in seq(x)) {
        dfm[i, ] <- apply(x[[i]], 2, mean, na.rm = TRUE)
        dfs[i, ] <- apply(x[[i]], 2, sd, na.rm = TRUE)
        x1 <- x[[i]][, 1]
        x2 <- x[[i]][, 2]
        pvalues[i] <- t.test(x1[!is.na(x1)], x2[!is.na(x2)])$p.value
    }
    
    # change data structure
    dfm <- t(dfm)
    dfs <- t(dfs)
    
    # calculate the y-coordinates for plotting *
    maxy <- max(dfm + dfs) * 1.1
    stepy <- max(dfm + dfs) * 0.1

    # bar chart
    bb <- barplot(dfm, beside = TRUE, ylim = c(0, maxy + 2 * stepy), col = col, border = col, xlab = xlab, ylab = ylab)
    
    # error bar
    arrows(bb, dfm - dfs, bb, dfm + dfs, code = 3, lwd = 1, angle = 90, length = 0.25 / length(sample.labels))
     
    # write *
    for (i in 1:length(pvalues)) {
        xi <- bb[, i]
        yi <- dfm[, i] + stepy * 1.5
        maxyi <- max(yi) + stepy
        if (pvalues[i] < 0.05) {
            lines(c(xi[1], xi[1], xi[2], xi[2]), c(yi[1], maxyi, maxyi, yi[2]))
            if (pvalues[i] < 0.01) {
                text((xi[1] + xi[2]) / 2, maxyi + stepy / 4, "**")
            } else if (pvalues[i] < 0.05) {
                text((xi[1] + xi[2]) / 2, maxyi + stepy / 4, "*")
            }
        }
    }
    
    # graph legend
    legend("topright", legend = condition.labels, fill = col, col = col, border = col, box.lwd = 0, box.lty = 0)
}

#####

ggpairs_2 = function(data = iris[,-ncol(iris)], group = iris$Species,title="TITLE", PDF_file = "tmp.pdf", pdf_h=5, pdf_w=5){
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(GGally))
    
    N_col = ncol(data)
    ggp   = ggpairs(data, upper='blank', diag='blank', lower='blank',title=title)
    
    for(iii in 1:N_col) {
      x = data[,iii]
      tmp_df = data.frame(x, gr = group)
      tmp_df = tmp_df[!apply(tmp_df,1,anyNA),]
      x=tmp_df$x
      p = ggplot(tmp_df, aes(x)) +
           theme_bw()+
           theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
      
      if (class(x) == 'factor') {
        p = p + geom_bar(aes(fill=gr), color='grey20')
      } else {
        bw = (max(x)-min(x))/10
        p  = p + geom_histogram(binwidth=bw, aes(fill=gr), color='grey20')+
              geom_line(eval(bquote(aes(y=..count..*.(bw)))), stat='density')
      }
    
      p   = p + geom_label(data=data.frame(x=-Inf, y=Inf, label=colnames(data)[iii]), aes(x=x, y=y, label=label), hjust=0, vjust=1)
      ggp = putPlot(ggp, p, iii, iii)
    }
    
    zcolat <- seq(-1, 1, length=81)
    zcolre <- c(zcolat[1:40]+1, rev(zcolat[41:81]))
    
    for(iii in 1:(N_col-1)) {
      for(jjj in (iii+1):N_col) {
        x = as.numeric(data[,iii])
        y = as.numeric(data[,jjj])
        r = cor(x, y, method='spearman', use='pairwise.complete.obs') %>% as.numeric_f()
        zcol = lattice::level.colors(r, at=zcolat,
                                     col.regions=colorRampPalette(c(scales::muted('red'), 'white', scales::muted('blue')), space='rgb')(81))
        textcol = ifelse(abs(r) < 0.4, 'grey20', 'white')
        ell = ellipse::ellipse(r, level=0.95, type='l', npoints=50, scale=c(.2, .2), centre=c(.5, .5))
        p   = ggplot(data.frame(ell), aes(x=x, y=y)) +
               theme_bw() + theme(plot.background=element_blank(),
                                  panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                  panel.border=element_blank(), axis.ticks=element_blank()           ) + 
               geom_polygon(fill=zcol, color=zcol) + 
               geom_text(data=NULL, x=.5, y=.5, label=100*round(r, 2), size=6, col=textcol)
        ggp = putPlot(ggp, p, iii, jjj)
      }
    }

    for(jjj in 1:(N_col-1)) {
      for(iii in (jjj+1):N_col) {
        x = data[,jjj]
        y = data[,iii]
        tmp_df = data.frame(x, y, gr=group)
        tmp_df = tmp_df[!apply(tmp_df,1,anyNA),]
        p = ggplot(tmp_df, aes(x=x, y=y, color=gr))+
             theme_bw() + 
             theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
        if (class(x) == 'factor') {
          p = p + geom_boxplot(aes(group=x), alpha=3/6, outlier.size=0, fill='white') +
                  geom_point(position=position_jitter(w=0.4, h=0), size=1)
        } else {
          p = p + geom_point(size=1)
        }
        ggp = putPlot(ggp, p, iii, jjj)
      }
    }
    
    pdf_3(PDF_file,h=pdf_h,w=pdf_w)
     print(ggp)
    dev.off()
}

####

ggpairs_3 = function(data = iris[,-ncol(iris)], group = iris$Species,title="TITLE", PDF_file = "tmp.pdf", pdf_h=5, pdf_w=5){
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(GGally))
    
    N_col = ncol(data)
    ggp   = ggpairs(data, upper='blank', diag='blank', lower='blank',title=paste0(today()," ",title))
    
    for(iii in 1:N_col) {
      x = data[,iii]
      tmp_df = data.frame(x, gr = group)
      tmp_df = tmp_df[!apply(tmp_df,1,anyNA),]
      x=tmp_df$x
      p = ggplot(tmp_df, aes(x)) +
           theme_bw()+
           theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
      
      if (class(x) == 'factor') {
        p = p + geom_bar(aes(fill=gr), color='grey20')
      } else {
        bw = (max(x)-min(x))/10
        p  = p + geom_histogram(binwidth=bw, aes(fill=gr), color='grey20')+
              geom_line(eval(bquote(aes(y=..count..*.(bw)))), stat='density')
      }
    
      p   = p + geom_label(data=data.frame(x=-Inf, y=Inf, label=colnames(data)[iii]), aes(x=x, y=y, label=label), hjust=0, vjust=1)
      ggp = putPlot(ggp, p, iii, iii)
    }
    
    zcolat <- seq(-1, 1, length=81)
    zcolre <- c(zcolat[1:40]+1, rev(zcolat[41:81]))
    
    for(iii in 1:(N_col-1)) {
      for(jjj in (iii+1):N_col) {
        x = as.numeric(data[,iii])
        y = as.numeric(data[,jjj])
        r = cor(x, y, method='spearman', use='pairwise.complete.obs') %>% as.numeric_f()
        zcol = lattice::level.colors(r, at=zcolat,
                                     col.regions=colorRampPalette(c(scales::muted('red'), 'white', scales::muted('blue')), space='rgb')(81))
        textcol = ifelse(abs(r) < 0.4, 'white', 'red')
        ell = ellipse::ellipse(r, level=0.95, type='l', npoints=50, scale=c(.2, .2), centre=c(.5, .5))
        p   = ggplot(data.frame(ell), aes(x=x, y=y)) +
               theme_bw() + theme(plot.background=element_blank(),
                                  panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                  panel.border=element_blank(), axis.ticks=element_blank()           ) + 
               geom_polygon(fill=zcol, color=zcol) + 
               geom_text(data=NULL, x=.5, y=.5, label=100*round(r, 2), size=6, col=textcol)
        ggp = putPlot(ggp, p, iii, jjj)
      }
    }

    for(jjj in 1:(N_col-1)) {
      for(iii in (jjj+1):N_col) {
        x = data[,jjj]
        y = data[,iii]
        tmp_df = data.frame(x, y, gr=group)
        tmp_df = tmp_df[!apply(tmp_df,1,anyNA),]
        p = ggplot(tmp_df, aes(x=x, y=y, color=gr))+
             theme_bw() + 
             theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
        if (class(x) == 'factor') {
          p = p + geom_boxplot(aes(group=x), alpha=3/6, outlier.size=0, fill='white') +
                  geom_point(position=position_jitter(w=0.4, h=0), size=1)
        } else {
          p = p + geom_point(size=1)
        }
        ggp = putPlot(ggp, p, iii, jjj)
      }
    }

    if(PDF_file=="tmp.pdf"){PDF_file=paste0(today(),"_",title,".pdf")}
    
    pdf_3(PDF_file,h=pdf_h,w=pdf_w)
     print(ggp)
    dev.off()
}

####

Greek_letter_convert = function(original){gsub("α","a",original) %>%
                                           gsub("β","b",.) %>%
                                           gsub("γ","c",.) %>%
                                           gsub("κ","k",.) %>%
                                           gsub("θ","q",.) %>%
                                           gsub("σ","s",.)}

####

plot_theme    = function(size=10,theme="classic",legend=TRUE,x_angle=90){
                                     if(theme == "classic"){p = theme_classic()}
                                     if(theme == "void"){   p = theme_void()}
                                     if(x_angle==90){p =  p + theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))}else{
                                                     p =  p + theme(axis.text.x  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = x_angle , hjust = 1 , vjust = 1 ))
                                                    }
                                     p = p + theme(axis.text.y  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
                                                   axis.title.x = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
                                                   axis.title.y = element_text(colour = "black", size = size, face = "bold", family = "Helvetica", angle = 90),
                                                   plot.title   = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))
                                     if(legend){ p = p + theme(legend.title = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"),
                                                               legend.text  = element_text(colour = "black", size = size, face = "bold", family = "Helvetica"))}else{
                                                 p = p + theme(legend.position  = "none")
                                               }
                                     return(p)
                                    }

####

Volcano_DEG_label = function(data, celltype, label_up = NULL, label_down = NULL, log2FC = 0.5, thresh = 1e-300){
  data1 = dplyr::filter(data, abs(avg_log2FC) > log2FC) %>% dplyr::filter(., p_val_adj < 0.05)
  data2 = data[setdiff(rownames(data),rownames(data1)),]
  data1$label = "good"
  data2$label = "bad"
  data = rbind(data1,data2)
  
  data$p_val_adj[(data$p_val_adj<thresh)&(data$p_val_adj>=0)] = thresh
  
  data$log10pVal = -log10(data$p_val_adj)
  data$DEG = "ns"
  data$DEG[(data$avg_log2FC>log2FC)&(data$p_val_adj<0.05)&(data$label=="good")] = "up"
  data$DEG[(data$avg_log2FC< -log2FC)&(data$p_val_adj<0.05)&(data$label=="good")] = "down"
  data$DEG = data$DEG %>% factor(.,levels=c("up","down","ns"))
  
  data3 = dplyr::filter(data, DEG %in% c("up","down"))
  data3$label = rownames(data3)
  data4 = dplyr::filter(data, DEG %in% c("ns"))
  data4$label = NA
  data = rbind(data3,data4)

  up_genes <- data %>% dplyr::filter(., DEG %in% c("up"))
  down_genes <- data %>% dplyr::filter(., DEG %in% c("down"))
  sig_genes <- data %>% dplyr::filter(., DEG %in% c("up","down"))
  cols <- c("up" = "firebrick", "down" = "steelblue", "ns" = "grey") 
  
    if (!is.null(label_up) || !is.null(label_down)){
    labels = c(label_up, label_down)
    sig_genes$label2 = NA
    for(tmp_num in 1:nrow(sig_genes)){
      if(sig_genes$label[tmp_num] %in% labels){
          sig_genes$label2[tmp_num] = sig_genes$label[tmp_num]
    }
   }
  } else {
    sig_genes$label2 = sig_genes$label
  }

  up_genes_unlabel = sig_genes %>% dplyr::filter(DEG == "up" & is.na(label2))
  up_genes_label = sig_genes %>% dplyr::filter(DEG == "up" & !is.na(label2))
  down_genes_unlabel = sig_genes %>% dplyr::filter(DEG == "down" & is.na(label2))
  down_genes_label = sig_genes %>% dplyr::filter(DEG == "down" & !is.na(label2))

  # Calculate density near the avg_log2FC thresholds (e.g., ±0.5)
  density_threshold <- 0.1 # example density threshold for changing line transparency
  near_threshold_up <- sum(sig_genes$avg_log2FC > 0.4 & sig_genes$avg_log2FC < 0.6) / length(sig_genes$avg_log2FC)
  near_threshold_down <- sum(sig_genes$avg_log2FC > -0.6 & sig_genes$avg_log2FC < -0.4) / length(sig_genes$avg_log2FC)
  
  # Calculate alpha based on density (simplified example)
  alpha_vlines <- ifelse(near_threshold_up > density_threshold | near_threshold_down > density_threshold, 0.2, 0.8)

  p = ggplot(data,aes(x=avg_log2FC,y=log10pVal))+
   geom_point_rast(aes(colour = DEG), 
              alpha = 0.5, 
              shape = 16,
              size = 0.2,
              colour = 'grey') +
   geom_point_rast(data = up_genes_unlabel, 
              alpha = 0.5, 
              shape = 16,
              size = 0.5,
              colour = "firebrick") +
   geom_point_rast(data = down_genes_unlabel, 
              alpha = 0.5, 
              shape = 16,
              size = 0.5,
              colour = "steelblue") +
   geom_point_rast(data = up_genes_label,
              shape = 21,
              size = 1, 
              fill = "firebrick", 
              colour = "black") + 
   geom_point_rast(data = down_genes_label,
              shape = 21,
              size = 1, 
              fill = "steelblue", 
              colour = "black") + 
   geom_text_repel(data = sig_genes,aes(label = label2),size = 4)+ 
   theme_classic()+
   ylim(c(NA, 320)) +
   scale_x_continuous()+
   theme(axis.text.x=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
          axis.text.y=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
              axis.title=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
              plot.title=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
              legend.title=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
               legend.text=element_text(colour="black",size=13,face="bold", family = "Helvetica"))+
       labs(title=paste0(celltype),x="avg_log2FC",y="−log10(adjusted_pValue)",fill="DEG") +
   geom_vline(xintercept=log2FC,col="black",size=0.5,linetype="dashed", alpha=alpha_vlines)+
   geom_vline(xintercept=-log2FC,col="black",size=0.5,linetype="dashed", alpha=alpha_vlines)+
   geom_hline(yintercept=1.3,col="black",size=0.5,linetype="dashed", alpha=alpha_vlines)

   return(p)
}

####

Volcano_DEG_label_v2 = function(data, celltype, label_up = NULL, label_down = NULL, xlim = 4, ylim = 150, log2FC = 0.5, p_val_thresh = 0.05, thresh = 1e-150){
  data1 = dplyr::filter(data, abs(avg_log2FC) > log2FC) %>% dplyr::filter(., p_val_adj < p_val_thresh)
  data2 = data[setdiff(rownames(data),rownames(data1)),]
  data1$label = "good"
  data2$label = "bad"
  data = rbind(data1,data2)
  
  data$p_val_adj[(data$p_val_adj<thresh)&(data$p_val_adj>=0)] = thresh
  
  data$log10pVal = -log10(data$p_val_adj)
  data$DEG = "ns"
  data$DEG[(data$avg_log2FC>log2FC)&(data$p_val_adj<p_val_thresh)&(data$label=="good")] = "up"
  data$DEG[(data$avg_log2FC< -log2FC)&(data$p_val_adj<p_val_thresh)&(data$label=="good")] = "down"
  data$DEG = data$DEG %>% factor(.,levels=c("up","down","ns"))
  
  data3 = dplyr::filter(data, DEG %in% c("up","down"))
  data3$label = rownames(data3)
  data4 = dplyr::filter(data, DEG %in% c("ns"))
  data4$label = NA
  data = rbind(data3,data4)

  up_genes <- data %>% dplyr::filter(., DEG %in% c("up"))
  down_genes <- data %>% dplyr::filter(., DEG %in% c("down"))
  sig_genes <- data %>% dplyr::filter(., DEG %in% c("up","down"))
  cols <- c("up" = "firebrick", "down" = "steelblue", "ns" = "grey") 
  
    if (!is.null(label_up) || !is.null(label_down)){
    labels = c(label_up, label_down)
    sig_genes$label2 = NA
    for(tmp_num in 1:nrow(sig_genes)){
      if(sig_genes$label[tmp_num] %in% labels){
          sig_genes$label2[tmp_num] = sig_genes$label[tmp_num]
    }
   }
  } else {
    sig_genes$label2 = sig_genes$label
  }

  up_genes_unlabel = sig_genes %>% dplyr::filter(DEG == "up" & is.na(label2))
  up_genes_label = sig_genes %>% dplyr::filter(DEG == "up" & !is.na(label2))
  down_genes_unlabel = sig_genes %>% dplyr::filter(DEG == "down" & is.na(label2))
  down_genes_label = sig_genes %>% dplyr::filter(DEG == "down" & !is.na(label2))

  # Calculate density near the avg_log2FC thresholds (e.g., ±0.5)
  density_threshold <- 0.1 # example density threshold for changing line transparency
  near_threshold_up <- sum(sig_genes$avg_log2FC > 0.4 & sig_genes$avg_log2FC < 0.6) / length(sig_genes$avg_log2FC)
  near_threshold_down <- sum(sig_genes$avg_log2FC > -0.6 & sig_genes$avg_log2FC < -0.4) / length(sig_genes$avg_log2FC)
  
  # Calculate alpha based on density (simplified example)
  alpha_vlines <- ifelse(near_threshold_up > density_threshold | near_threshold_down > density_threshold, 0.2, 0.8)

  p = ggplot(data,aes(x=avg_log2FC,y=log10pVal))+
   geom_point_rast(aes(colour = DEG), 
              alpha = 0.5, 
              shape = 16,
              size = 0.2,
              colour = 'grey') +
   geom_point_rast(data = up_genes_unlabel, 
              alpha = 0.5, 
              shape = 16,
              size = 0.5,
              colour = "firebrick") +
   geom_point_rast(data = down_genes_unlabel, 
              alpha = 0.5, 
              shape = 16,
              size = 0.5,
              colour = "steelblue") +
   geom_point_rast(data = up_genes_label,
              shape = 21,
              size = 1, 
              fill = "firebrick", 
              colour = "black") + 
   geom_point_rast(data = down_genes_label,
              shape = 21,
              size = 1, 
              fill = "steelblue", 
              colour = "black") + 
   geom_text_repel(data = sig_genes,aes(label = label2),size = 2)+ 
   theme_classic()+
   ylim(c(NA, ylim)) + xlim(c(-xlim,xlim)) +
   scale_x_continuous()+
   theme(axis.text.x=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
          axis.text.y=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
              axis.title=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
              plot.title=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
              legend.title=element_text(colour="black",size=13,face="bold", family = "Helvetica"),
               legend.text=element_text(colour="black",size=13,face="bold", family = "Helvetica"))+
       labs(title=paste0(celltype),x="avg_log2FC",y="−log10(adjusted_pValue)",fill="DEG") +
   geom_vline(xintercept=log2FC,col="black",size=0.5,linetype="dashed", alpha=alpha_vlines)+
   geom_vline(xintercept=-log2FC,col="black",size=0.5,linetype="dashed", alpha=alpha_vlines)+
   geom_hline(yintercept=-log10(p_val_thresh),col="black",size=0.5,linetype="dashed", alpha=alpha_vlines)

   return(p)
}

####
Get_Heatmap = function(df, genes_list, colors){
  ht <- Heatmap(
  t(df[genes_list,]),                  
  name = "Expression",            
  cluster_rows = FALSE,           
  cluster_columns = TRUE,         
  show_row_names = TRUE,          
  show_column_names = TRUE,       
  col = colors,                  
  row_names_side = "left",        
  column_names_side = "top",      
  heatmap_legend_param = list(
    title = "Expression",
    at = c(-2, 0, 2),             
    labels = c("-2", "0", "2")
    )
  )

  return(ht)
}


####
plot_significances <- function(
    significances,
    p_key,
    value_key,
    group_key,
    enrichment_key = "enrichment",
    enriched_label = "enriched",
    pmax = 0.05,
    pmin = 1e-5,
    annotate_pvalues = TRUE,
    annotation_fontsize = 8,
    value_cluster = FALSE,
    group_cluster = FALSE,
    value_order = NULL,
    group_order = NULL,
    method = method
) { library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(reshape2)
    
    # Constants
    small_value <- 1e-300
    max_log <- -log10(pmin)
    min_log <- -log10(pmax)
    
    # Validate required columns in significances
    required_columns <- c(p_key, value_key, group_key)
    missing_columns <- setdiff(required_columns, colnames(significances))
    if (length(missing_columns) > 0) {
        stop("The following required columns are missing in 'significances': ", 
             paste(missing_columns, collapse = ", "), ".")
    }
    
    # Validate enrichment key
    depleted_label <- NULL
    if (!is.null(enrichment_key)) {
        if (!enrichment_key %in% colnames(significances)) {
            stop(paste0(
                'The column "', enrichment_key, '" does not exist in the supplied dataframe. ',
                'Set enrichment_key = NULL if this is intentional.'
            ))
        }
        
        unique_significance_labels <- unique(significances[[enrichment_key]])
        if (length(unique_significance_labels) == 1) {
            enriched_label <- unique_significance_labels[1]
        } else if (length(unique_significance_labels) > 2 || !(enriched_label %in% unique_significance_labels)) {
            stop(paste0(
                'The column "', enrichment_key, 
                '" must have exactly 2 values: "', enriched_label, 
                '" and another (e.g., "depleted"). Found: ', 
                paste(unique_significance_labels, collapse = ", "), '.'
            ))
        } else {
            depleted_label <- setdiff(unique_significance_labels, enriched_label)[1]
        }
    }
    
    # Prepare matrices for plotting
    if (!is.null(depleted_label)) {
        enr_e <- reshape2::acast(significances[significances[[enrichment_key]] == enriched_label, ], 
                                 formula = paste(value_key, group_key, sep = "~"), 
                                 value.var = p_key)
        enr_p <- reshape2::acast(significances[significances[[enrichment_key]] == depleted_label, ], 
                                 formula = paste(value_key, group_key, sep = "~"), 
                                 value.var = p_key)
        
        enr_e[is.na(enr_e)] <- small_value
        enr_p[is.na(enr_p)] <- small_value
        
        enr <- ifelse(enr_e < enr_p, -log10(enr_e), log10(enr_p))
        ann <- ifelse(enr_e < enr_p, enr_e, enr_p)
    } else {
        ann <- reshape2::acast(significances, 
                               formula = paste(value_key, group_key, sep = "~"), 
                               value.var = p_key)
        enr <- -log10(ann)
    }
    
    # Filter insignificant values
    enr[ann > pmax] <- 0
    ann[ann > pmax] <- NA
    
    # Reorder rows and columns if specified
    if (!is.null(value_order)) {
        if (anyNA(value_order)) stop("value_order contains NA values. Please remove or replace them.")
        enr <- enr[value_order, , drop = FALSE]
    }
    if (!is.null(group_order)) {
        if (anyNA(group_order)) stop("group_order contains NA values. Please remove or replace them.")
        enr <- enr[, group_order, drop = FALSE]
    }
    
    # Define color scale
    depleted_color1 <- rgb(0.30196078431372547, 0.5215686274509804, 0.7098039215686275)
    depleted_color2 = rgb(0.7803922, 0.8243137, 0.8619608)
    null_color <- rgb(0.9,0.9,0.9)
    enriched_color1 <- rgb(1.0, 0.07058823529411765, 0.09019607843137255)
    enriched_color2 <- rgb(0.92, 0.7341176, 0.7380392)

    # Custom color function
    custom_color_fun <- function(values) {
        sapply(values, function(value) {
            if (value > -min_log & value < min_log) {
                return(null_color)
            }
            colorRamp2(
                breaks = c(-max_log, -min_log, min_log, max_log),
                colors = c(depleted_color1, depleted_color2, enriched_color2, enriched_color1)
            )(value)
        })
    }
    
    # Replace col_fun with the custom logic
    col_fun <- custom_color_fun

    # Configure clustering with "average" method
    row_dend <- if (value_cluster) hclust(dist(enr), method = method) else FALSE
    col_dend <- if (group_cluster) hclust(dist(t(enr)), method = method) else FALSE
    
    # Create heatmap
    heatmap <- ComplexHeatmap::Heatmap(
        matrix = enr,
        name = "-log10(p-value)",
        col = col_fun,
        cluster_rows = row_dend,
        cluster_columns = col_dend,
        show_row_dend = !is.logical(row_dend),
        show_column_dend = !is.logical(col_dend),
        heatmap_legend_param = list(
            at = c(-min_log, min_log, -max_log, max_log),
            labels = c(-pmax, pmax, -pmin, pmin)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (annotate_pvalues && !is.na(ann[i, j])) {
                grid.text(sprintf("%.2e", ann[i, j]), x, y, gp = gpar(fontsize = annotation_fontsize))
            }
        }
    )
  
    # Draw heatmap
    draw(heatmap)
}

####
rank_norm <- function(df){  ## df: celltype*Pts
  rank_norm_ = function(x) {
  ranks <- rank(x, ties.method = "average")  # Rank the values in each column
  normalized_values <- qnorm((ranks - 0.5) / length(ranks))  # Apply the inverse normal transformation
  return(normalized_values)
  }
  df_norm <- as.data.frame(apply(df, 1, rank_norm_)) %>% t()
  return(df_norm)
}

####
clr_norm = function(df){        # df: celltype*Pts
  df[df == 0] <- 1e-6
  df_clr <- as.data.frame(apply(df, 1, function(x) clr(x))) %>% t()
  return(df_clr)
}

####
plot_norm_hist = function(df_orig, df_norm, path){  # df: celltype*Pts
  p = list()
  for(tmp_num in 1:nrow(df_orig)){
      df_tmp = data.frame( original  = df_orig[tmp_num,] %>% as.numeric(),
                           normalize = df_norm[tmp_num,] %>% as.numeric(),
                           name      = rownames(df_orig)[tmp_num])
      p_1 =  ggplot(df_tmp , aes_string(x = "original", fill = "name")) +
                geom_density(alpha = .3, col = "black") +
                scale_fill_manual(values=c("red")) + 
                plot_theme()+
                theme(legend.position="none")+
                theme(axis.text.y=element_blank())+
                labs(title=rownames(df_orig)[tmp_num])
      p_2 =  ggplot(df_tmp , aes_string(x = "normalize", fill = "name")) +
                geom_density(alpha = .3, col = "black") +
                scale_fill_manual(values=c("blue")) + 
                plot_theme()+
                theme(legend.position="none")+
                theme(axis.text.y=element_blank())
      p[[tmp_num]] = p_1 + p_2 + plot_layout(ncol=1)+
           plot_annotation(paste0("UC_subset_proportion_correlation"," ",rownames(df_orig)[tmp_num]),
                           theme = theme(plot.title = element_text(size = 15, hjust = 0.5,face="bold", family = "Helvetica")))
  }
  p_sum = patchwork::wrap_plots(p, ncol = 8)
  pdf_3(path,h=12,w=16)
    plot(p_sum)
  dev.off()
}

####
plot_corrmatrix = function(df, col_list,
                           cor_method = 'spearman',
                           clustering_method = 'ward.D2',
                           save_path){
  df_cor = cor(df %>% as.matrix() %>% t() ,method=cor_method)
  col = colorRampPalette(col_list)
  pheatmap_result = pheatmap::pheatmap(df_cor,
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE , 
                   cellwidth = 10, cellheight = 10,
                   color = col(50),
                   clustering_method = clustering_method,
                   breaks = seq(-1,1,length.out=50),
                   filename=save_path,
                   border_color = "black" ,na_col = "gray90")
  return(pheatmap_result)
}


####
TACCO_neighbour_enrichment_boxplot <- function(group_A, group_B, data_dir, colors){
    samples <- c(group_A, group_B)
    dot_fill <- colors
    box_fill <- colors
    
    # Read and format data
    z_values_all <- map_dfr(samples, function(s) {
      df <- read.delim(file.path(data_dir, paste0(s, "_z.txt")), row.names = 1, check.names = FALSE)
      df %>%
        rownames_to_column("from") %>%
        dplyr::filter(from == "LUM") %>%
        pivot_longer(-from, names_to = "celltype", values_to = "z") %>%
        mutate(
          sample = s,
          group = ifelse(s %in% group_A, "Group A", "Group B")
        )
    })
    
    # Calculate stats per celltype
    stats_df <- z_values_all %>%
      group_by(celltype) %>%
      summarise(
        pval = wilcox.test(z ~ group, exact = FALSE)$p.value,
        y_line = max(z, na.rm = TRUE) + diff(range(z, na.rm = TRUE)) * 0.1
      ) %>%
      mutate(pval_label = paste0("p = ", signif(pval, 3)))
    
    # Merge back
    z_plot <- z_values_all %>%
      left_join(stats_df, by = "celltype")
    
    # Plot
    p <- ggplot(z_plot, aes(x = group, y = z)) +
      geom_boxplot(aes(fill = group), width = 0.3, color = "black", alpha = 0.4) +
      geom_quasirandom(
        aes(fill = group),
        shape = 21, size = 5, stroke = 1, color = "black",
        width = 0.15, dodge.width = 0.5, alpha = 0.8
      ) +
      geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line), inherit.aes = FALSE) +
      geom_segment(aes(x = 1, xend = 1, y = y_line - 0.05, yend = y_line), inherit.aes = FALSE) +
      geom_segment(aes(x = 2, xend = 2, y = y_line - 0.05, yend = y_line), inherit.aes = FALSE) +
      geom_text(aes(x = 1.5, y = y_line + 0.05 * max(y_line + 1, na.rm = TRUE), label = pval_label),
                inherit.aes = FALSE, size = 4) +
      scale_fill_manual(values = box_fill) +
      facet_wrap(~celltype, scales = "free_y") +
      labs(y = "Z-score", x = NULL) +
      theme_classic(base_size = 14) +
      theme(
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none"
      )

      return (p)
}

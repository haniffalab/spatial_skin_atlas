library(Seurat)
library(visutils)
setwd("/nfs/cellgeni/pasham/projects/2211.adult.skin/figures/paper_figures/")

# load data ############
gene.descr=readRDS('data/gene.descr.rds')
# list of (for different annotations) distance to dermis/epidermis junction * celltype * sample summary matrices
# for now there is just one annotation: lvl4
dfsmtx.c2l = readRDS('data/dfsmtx.c2l.raw.rds')

# distance to dermis/epidermis junction * gene * sample summary matrix (lcpm)
dfsmtx.ge = readRDS('data/dfsmtx.ge.all.rds')

# sample metadata
meta = readRDS('data/meta_filtered.rds')
all(rownames(meta) == dimnames(dfsmtx.c2l)[[3]])
all(rownames(meta) == dimnames(dfsmtx.ge)[[3]])

# make sample groups to be used in paper
meta$groups01 = '-'
meta$groups01[meta$State == 'EL'] = 'AD L'
meta$groups01[meta$dataset_state == 'bayanne_PL'] = 'PSO L'
meta$groups01[meta$dataset_state == 'bayanne_PnL'] = 'PSO NL'
meta$groups01[meta$dataset_state == 'bayanne_EnL'] = 'AD NL'
meta$groups01[meta$dataset_state == 'bayanne_H'] = 'H'

table(meta$groups01,meta$dataset_state)

# celltype info
celltypes = read.csv('data/celltypes.csv',row.names = 1)
celltypes = celltypes[nrow(celltypes):1,]
lvl0_cols = char2col(celltypes$lvl0)



# cond_cols = c('AD NL'='yellow',
#               'PSO NL'='#fdcac9',
#               'AD L'='#FF7F00',
#               'PSO L'='darkred')

cond_cols = c('AD NL'='#CCCAEA',
              'PSO NL'='#F2CECE',
              'AD L'='#8D8BBC',
              'PSO L'='#DD9B9B')

path2figs = '.'
dir.create(path2figs,recursive = TRUE)
sd.mult = 2

# rename celltypes
ct_rename = c('KC_HF: IRS_cycling' = 'KC_HF: Matrix',
              'KC_HF: ORS-Bulb_infundibulum' = 'KC_HF: SPON2+')
f = celltypes$lvl4_annotation %in% names(ct_rename)
celltypes[f,'lvl4_annotation'] = ct_rename[celltypes[f,'lvl4_annotation']]

f = dimnames(dfsmtx.c2l$lvl4)[[2]] %in% names(ct_rename)
dimnames(dfsmtx.c2l$lvl4)[[2]][f] = ct_rename[dimnames(dfsmtx.c2l$lvl4)[[2]][f]]

# reorder celltypes
# KS order: 
kc_order = c('KC1', 'KC1-2_cycling', 'KC2', 'KC3', 'KC2_3_cycling', 'KC4', 'KC5', 'KCinflamm_basal', 'KCinflamm_int_late')
# other alphabetically
f = which(celltypes$lvl0 != 'lloyd.man.ann')
lvl0_depth = tapply(celltypes$mean_depth, celltypes$lvl0, mean)
celltypes[f,] = celltypes[f[order(-lvl0_depth[celltypes$lvl0[f]],celltypes$lvl4_annotation[f])],]
rownames(celltypes) = celltypes$lvl4_annotation

# reorder KC
f = which(celltypes$lvl0 == 'KC')
celltypes[f,] = celltypes[kc_order,]

all(rownames(celltypes) == celltypes$lvl4_annotation)
all(dimnames(dfsmtx.c2l$lvl4)[[2]] %in% celltypes$lvl4_annotation)

# plot figures #########

# _01 fig3_visutils_validation (3d) ###############
# distance range to use
dist_range = as.character(-15:5)
celltypes2plot = c('Treg','T_Prolif','pDC','Treg_LRRC32+')
ann = 'lvl4'

cairo_pdf(paste0(path2figs,'/01_fig3_visutils_validation_3d.pdf'),w=length(celltypes2plot)*2.2+2.2,h=2.2,onefile = TRUE)
par(mfrow=c(1,length(celltypes2plot)+1),mar=c(2,3,1,0),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,0,0))
for(ct in celltypes2plot){
  par(xaxt='n',yaxt='n')
  plotConditionsProfiles(dfsmtx.c2l[[ann]][dist_range,,],feature=ct,meta$groups01,
                         cols = cond_cols,lwd=3,sd.mult = sd.mult,main=ct,
                         legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
  par(xaxt='l')
  axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
  abline(v=0,lty=3)
}
plot.new()
legend('topleft',bty = 'n',lwd=2,col=cond_cols,legend = names(cond_cols))
dev.off()

# _02 supp_figure_visutils_locations_steadystate (S13) ############
table(meta$State,meta$dataset)
groups = c('AD NL' = 'Eczema non-lesional',
           'PSO NL' = 'Psoriasis non-lesional',
           'H'   = 'Healthy (no inflammatory skin disease')
ann = 'lvl4'
# distance range to use
dist_range = as.character(-25:3)

mtxs = lapply(names(groups),function(s){
  t(apply(dfsmtx.c2l[[ann]][dist_range,,meta$groups01==s],1:2,mean,na.rm=TRUE))
})
names(mtxs) = names(groups)

#cmn_maxs = apply(sapply(mtxs,function(d)apply(d,1,max,na.rm=T)),1,max)
cairo_pdf(paste0(path2figs,'/02_supp_figure_visutils_locations_steadystate_S13.pdf'),w=10,h=13,onefile = TRUE)
par(mfrow=c(3,1),bty='n',mar=c(14,1,1,4))
for(n in names(groups)){
  d = mtxs[[n]]
  # per condition maxes
  d = sweep(d,1,apply(d,1,max,na.rm=T),'/')
  # common maxes
  # d = sweep(d,1,cmn_maxs,'/')
  imageWithText(d[celltypes$lvl4_annotation,],'',yaxt='n',yaxlab=NULL,cex.axis.x = 0.8,
                colAnns = list(lvl0=celltypes$lvl0),
                colAnnCols = list(lvl0=lvl0_cols),
                main=groups[n])
  abline(h=which(colnames(d)=='0'),lty=2)
  axis(4,1:ncol(d),labels = colnames(d),las=2)
  mtext('Distance to epidermis−dermis interface (spots)',4,line=2.6,cex=0.7)
}
dev.off()

# _03 supp_figure_visutils_location_disease (S14) ############
table(meta$State,meta$dataset)
groups = c('AD L' = 'Eczema lesional combined',
           'PSO L' = 'Psoriasis lesional')
ann = 'lvl4'
# distance range to use
dist_range = as.character(-25:3)

mtxs = lapply(names(groups),function(s){
  t(apply(dfsmtx.c2l[[ann]][dist_range,,meta$groups01==s],1:2,mean,na.rm=TRUE))
})
names(mtxs) = names(groups)

#cmn_maxs = apply(sapply(mtxs,function(d)apply(d,1,max,na.rm=T)),1,max)
cairo_pdf(paste0(path2figs,'/03_supp_figure_visutils_location_disease_S14.pdf'),w=10,h=13,onefile = TRUE)
par(mfrow=c(3,1),bty='n',mar=c(14,1,1,4))
for(n in names(groups)){
  d = mtxs[[n]]
  # per condition maxes
  d = sweep(d,1,apply(d,1,max,na.rm=T),'/')
  # common maxes
  # d = sweep(d,1,cmn_maxs,'/')
  imageWithText(d[celltypes$lvl4_annotation,],'',yaxt='n',yaxlab=NULL,cex.axis.x = 0.8,
                colAnns = list(lvl0=celltypes$lvl0),
                colAnnCols = list(lvl0=lvl0_cols),
                main=groups[n])
  abline(h=which(colnames(d)=='0'),lty=2)
  axis(4,1:ncol(d),labels = colnames(d),las=2)
  mtext('Distance to epidermis−dermis interface (spots)',4,line=2.6,cex=0.7)
}
dev.off()

# _04 ext_figure_visutils_validation_extra (E5b) ###################
dist_range = as.character(-15:5)
celltypes2plot = c('Pericyte1','cDC2: MMP12hi','Tc2','Tc3_IFNGhi')
ann = 'lvl4'

setdiff(celltypes2plot,dimnames(dfsmtx.c2l$lvl4)[[2]])

cairo_pdf(paste0(path2figs,'/04_ext_figure_visutils_validation_extra_E5b.pdf'),w=length(celltypes2plot)*2.2+2.2,h=2.2,onefile = TRUE)
par(mfrow=c(1,length(celltypes2plot)+1),mar=c(2,3,1,0),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,0,0))
for(ct in celltypes2plot){
  par(xaxt='n',yaxt='n')
  plotConditionsProfiles(dfsmtx.c2l[[ann]][dist_range,,],feature=ct,meta$groups01,
                         cols = cond_cols,lwd=3,sd.mult = sd.mult,main=ct,
                         legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
  par(xaxt='l')
  axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
  abline(v=0,lty=3)
}
plot.new()
legend('topleft',bty = 'n',lwd=2,col=cond_cols,legend = names(cond_cols))
dev.off()

# _05 ext_figure_ccl1_and_ilc3 (E1f/g) ###################
dist_range = as.character(-15:5)

celltype = 'ILC3_CCL1+PTGDS+'
gene_name='CCL1'
gene_id = gene.descr$ens_id[gene.descr$name==gene_name]
ann = 'lvl4'


cairo_pdf(paste0(path2figs,'/05_ext_figure_ccl1_and_ilc3_E1fg.pdf'),w=2*2.2+2.2,h=2.2,onefile = TRUE)

par(mfrow=c(1,3),mar=c(2,3,1,0),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,0,0))
par(xaxt='n',yaxt='n')
plotConditionsProfiles(dfsmtx.ge[dist_range,,],feature=gene_id,meta$groups01,
                         cols = cond_cols,lwd=3,sd.mult = sd.mult,main=gene_name,
                         legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
par(xaxt='l')
axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
abline(v=0,lty=3)

plot.new()
legend('topleft',bty = 'n',lwd=2,col=cond_cols,legend = names(cond_cols))

par(xaxt='n')
plotConditionsProfiles(dfsmtx.c2l[[ann]][dist_range,,],feature=celltype,meta$groups01,
                       cols = cond_cols,lwd=3,sd.mult = sd.mult,main=celltype,
                       legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
par(xaxt='l')

axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
abline(v=0,lty=3)

dev.off()


# _06 ext_figure_eotaxins_and_neutrophil (E7a) ###################
dist_range = as.character(-15:5)

gene_names=c('CCL11','CCL24','CCL26','CXCL1','CXCL5','CXCL8')
genes = gene.descr[match(gene_names,gene.descr$name),]

# line is 0.2 in
cairo_pdf(paste0(path2figs,'/06_ext_figure_eotaxins_and_neutrophil_E7a.pdf'),w=3*2+1.6,h=2*2,onefile = TRUE)

par(mfrow=c(2,3),mar=c(2,3,3,0),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,0,8))
for(gene_id in genes$ens_id){
  par(xaxt='n',yaxt='n')
  plotConditionsProfiles(dfsmtx.ge[dist_range,,],feature=gene_id,meta$groups01,
                         cols = cond_cols,lwd=3,sd.mult = sd.mult,main=genes[gene_id,'name'],
                         legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
  par(xaxt='l')
  axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
  abline(v=0,lty=3)
}

legend(grconvertX(1,'npc','user'),grconvertY(1,'npc','user'),bty = 'n',lwd=2,col=cond_cols,legend = names(cond_cols),xpd=NA)

dev.off()

# _07 ext_figure_panel_for_stroma_stimulation (E7e) ###################
dist_range = as.character(-15:5)

gene_names=c('IL1B','TNF','IL17A','IL13')
genes = gene.descr[match(gene_names,gene.descr$name),]

# line is 0.2 in
cairo_pdf(paste0(path2figs,'/07_ext_figure_panel_for_stroma_stimulation_E7e.pdf'),w=2*2,h=2*2+1,onefile = TRUE)

par(mfrow=c(2,2),mar=c(2,3,3,0),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(5,0,0,0))
for(gene_id in genes$ens_id){
  par(xaxt='n',yaxt='n')
  plotConditionsProfiles(dfsmtx.ge[dist_range,,],feature=gene_id,meta$groups01,
                         cols = cond_cols,lwd=3,sd.mult = sd.mult,main=genes[gene_id,'name'],
                         legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
  par(xaxt='l')
  axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
  abline(v=0,lty=3)
}

legend(grconvertX(0.5,'ndc','user'),grconvertY(0,'nfc','user'),bty = 'n',lwd=2,col=cond_cols,legend = names(cond_cols),xpd=NA,ncol=2,xjust = 0.5)

dev.off()

# _08 cDC2: EREG+CCR7+ and LC  ###############
# distance range to use
dist_range = as.character(-15:5)
celltypes2plot = c('LC','cDC2: EREG+CCR7+')
ann = 'lvl4'

setdiff(celltypes2plot,colnames(dfsmtx.c2l[[ann]]))

cairo_pdf(paste0(path2figs,'/08_cDC2_EREG_CCR7_and_LC.pdf'),w=length(celltypes2plot)*2.2,h=2.2,onefile = TRUE)
par(mfrow=c(1,length(celltypes2plot)),mar=c(2,3,1,0),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,0,0),cex=0.66)
for(ct in celltypes2plot){
  par(xaxt='n',yaxt='n')
  plotConditionsProfiles(dfsmtx.c2l[[ann]][dist_range,,],feature=ct,meta$groups01,
                         cols = cond_cols,lwd=3,sd.mult = sd.mult,main=ct,
                         legend. = FALSE,xlim=c(-15,5),xlab='',ylab='',cilim = c(0,Inf))
  par(xaxt='l')
  axis(1,at=c(-10,0),lab=c('1 mm','d/epi border'),xpd=NA,tcl=-0.1,mgp=c(1.15,0.15,0))
  abline(v=0,lty=3)
}

legend('topleft',bty = 'n',lwd=2,col=cond_cols,legend = names(cond_cols))
dev.off()

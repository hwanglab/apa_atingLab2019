#!/usr/bin/env Rscript

#Authors: Chanjing Hong
#Date Created: 1/2/2018

library(data.table)
library(stringr)

jmessage <- function(msg2display,action='info',job_count=1) {
	now = Sys.time()
	if (action == 'error') {
		message(sprintf("[ERROR|%s] (%d) %s",now,job_count,msg2display))
		stop(msg2display)
	} else if (action == 'warning') {
		message(sprintf("[WARNING|%s] (%d) %s",now,job_count,msg2display))
	} else if (action == 'start' | action == 'banner') {
		message("=============================")
		message(sprintf("[TASK|%s] (%d) %s",now,job_count,msg2display))
	} else if (action == 'end') {
		message(sprintf("[TASK|%s] (%d) %s",now,job_count,msg2display))
		message("=============================")
	} else {
		message(sprintf("[INFO|%s] (%d) %s",now,job_count,msg2display))
	} 
}

print_table <- function(u_title,u_fontSize,ms3a) {
	
	title = textGrob(u_title,gp=gpar(fontsize=12))
	padding = unit(5,'mm')
	table = gtable_add_rows(tableGrob(ms3a,rows=NULL,theme=ttheme_default(base_size=u_fontSize)),
													heights = grobHeight(title) + padding,
													pos=0)
	
	table = gtable_add_grob(table,
													title,
													1,1,1,ncol(table),
													clip = 'off')
	grid.newpage()
	grid.draw(table)
}

file_filter_rows <- function(fn,colidx1,key) {
	fn_head = sprintf('%s.head',fn)
	system(sprintf('head -n1 %s > %s',fn,fn_head))
	fn_tmp = sprintf('%s.tmp',fn)
	system(sprintf('awk \'$%d == \"%s\" {print $0}\' %s > %s',colidx1,key,fn,fn_tmp))
	system(sprintf('cat %s >> %s',fn_tmp,fn_head))
	unlink(fn_tmp)
	fn_head
}

get_max_val_in_list <- function(my_list,absolute=F){
  M <- length(my_list)
  maxVals <- rep(-1e7,M)
  for (i in 1:M){
    maxVals[i] <- max(my_list[[i]],na.rm = T)
  }
  return(max(maxVals))
}


get_min_val_in_list <- function(my_list,absolute=F){
  M <- length(my_list)
  minVals <- rep(1e7,M)
  for (i in 1:M){
    minVals[i] <- min(my_list[[i]],na.rm = T)
  }
  return(min(minVals))
}

get_mount_dir <- function(){
  hostname <- Sys.getenv('HOSTNAME')
  if (hostname == 'lri-107577'){
    message('from author local directory/development purpose!')
    mount_prefix <- '/home/hongc2/orange_passport/apa_atingLab2019'
  } else {
    mount_prefix <- file.path(Sys.getenv('HOME'),'projects','apa_atingLab2019')
  }
  return(mount_prefix)
}

comma1k <- function(my_integers){
  
  my_integer_dc <- data.class(my_integers[1])
  
  if (my_integer_dc=='numeric'){
    char_integers_with_comma <- format(my_integers, big.mark=",", scientific=FALSE)
  } else if (my_integer_dc == 'character') {
    char_integers_with_comma <- format(as.numeric(my_integers), big.mark=",", scientific=FALSE)
  } else {
    stop(sprintf('cannot support the data type [%s]',my_integer_dc))
  }
  return(char_integers_with_comma)
}

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
# 
# b37_to_hg19_bam <- function(bam){
# 	chr_bam <- sprintf('%s.chr.bam',bam)
# 	cmd <- sprintf("samtools view -H %s | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | grep -vP '(SN:GL)|(SN:hs)' | samtools reheader - %s > %s",bam,bam,chr_bam)
# 	system(cmd)
# 	return(chr_bam)
# }


get_file_base <- function(filepath2) {
  message(filepath2)
  fbase <- basename(filepath2)
  if (endsWith(fbase,'.gz')) {
    fbase <- sub(".gz","",fbase)
  }
  
  #finally strip off the original file extension
  fbase0 <- tools::file_path_sans_ext(fbase)
  return(fbase0)
}

page_to_print <- function(range_str) {
  
  range2 <- strsplit(range_str,',')[[1]]
  range2a <- strsplit(range2,'-')
  range0 <- list()
  for (i in length(range2a)) {
    range3 <- range2a[[i]]
    range4 <- as.numeric(unlist(range3))
    if (length(range4)==2) {
      range0[[i]] <- range4[1]:range4[2]
    } else {
      range0[[i]] <- range4
    }
  }
  page_to_print <- unique(unlist(range0))
  return(page_to_print)
}
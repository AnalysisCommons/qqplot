LOCAL <- FALSE
if(LOCAL) {
  # data.files <- Sys.glob("~/DNAnexus/data/mb_test_data_chr*.Rda")
  data.files <- "~/DNAnexus/data/mb_test_data_chr22_color.Rda"
  p.col <- "p"
  filter <- "FALSE"
  use.lambda <- TRUE
  include.ci <- TRUE
  qq.color <- c("blue","red","yellow")
  qq.color.criteria <- c("maf<0.1","maf<0.05","maf<0.01")
  output.type <- "png"

  library(gap)

} else {

  ## Parse the arguments
  args <- commandArgs(trailingOnly=TRUE)
  args <- paste(args,collapse=" ")
  print(args)
  args <- unlist(strsplit(args,";"))
  print(args)

  ## Mandatory parameters
  # data.files <- args[7:length(args)] # old way
  data.files <- unlist(strsplit(args[8]," "))

  ## Optional parameters
  p.col <- args[1]
  filter <- args[2]
  use.lambda <- args[3]
  include.ci <- args[4]
  qq.color <- unlist(strsplit(args[5]," "))
  qq.color.criteria <- unlist(strsplit(args[6]," "))
  output.type <- args[7]

  ## Quick check
  if ((qq.color.criteria[1] != "FALSE") & (length(qq.color) != length(qq.color.criteria))) {
      stop("Must use the same number of colors as color filtering criteria. Colors: ",qq.color,", Criteria: ",qq.color.criteria)
  }
  if ((qq.color.criteria[1] != "FALSE") & (qq.color[1] == "color")) {
      stop("qq_color_criteria must be FALSE if you want to use your own custom 'color' data column.")
  }

  ## Install libraries
  install.packages("/gap_1.1-16.tar.gz", repos=NULL, type="source")
  library(gap)

}

## Functions
getLambda <- function(p) {
	chi <- qchisq(p,df=1,lower.tail=F)
	g <- median(chi,na.rm=TRUE)
	lambda <- round(g/0.456,digits=4)
	return(lambda)
}

qqtrunc <- function(p, threshold, ci, col, ...) {
        p <- na.omit(p)
        n <- length(p)
        p <- sort(p[p<=threshold])
        n.p <- length(p)

	if (ci)
	    c <- abs(qnorm(0.05/2))
        m <- (1:n)/(n + 1)
	xlims <- ylims <- rep(NA,2)
	xlims <- c(0, max(-log10(m[1:n.p]),na.rm=TRUE))
	ylims <- c(0, max(-log10(p),na.rm=TRUE))
        z <- qqplot(-log10(m[1:n.p]), -log10(p), xlab = "-log10(Expected)",
                ylab = "-log10(Observed)", col = qq.color, xlim=xlims, ylim=ylims, ...)
	if (ci) {
            v <- (1:n) * (n - (1:n) + 1)/(n + 1)^2/(n + 2)
            s <- sqrt(v)
            lcl <- m - c * s
            ucl <- m + c * s
            lid <- (lcl > 0)
            uid <- (ucl <= 1)
            a <- -log10(m[lid])
            b <- -log10(lcl[lid])
            c <- -log10(m[uid])
            d <- -log10(ucl[uid])
            
            lines(a, b)
            lines(c, d)
        }
        abline(0, 1, col = "red")
}
        
l.sing <- l.skat <- vector(mode="list",length=length(data.files))

## Added for custom colors
if (qq.color[1]=="color" | qq.color.criteria[1]!="FALSE") {
    l.col <- l.col.skat <- vector(mode="list",length=length(data.files))
} 

for (i in 1:length(data.files)) {
	cat("Loading ",data.files[i],"\n")
	load(data.files[i])
	
	# If using a filter, subset the data here
	
	# if (is.element(mac.col,names(sing))) {
	#	sing <- eval(parse(text = paste0("subset(sing,",mac.col,mac.filter,")")))
	#}

	if (filter != "FALSE") {
		message("Filtering singleSNP data")
		sing <- eval(parse(text = paste0("subset(sing,",filter,")")))
	}
	l.sing[[i]] <- sing[,p.col]
	l.skat[[i]] <- skat[,p.col]

        ## Added for custom colors
        if (qq.color[1]=="color") {
            l.col[[i]] <- sing[,"color"]
	    l.col.skat[[i]] <- skat[,"color"]
        } else if (qq.color.criteria[1]!="FALSE") {
            def.col <- ifelse(!is.element("black",qq.color), "black", ifelse(!is.element("blue",qq.color),"blue","gray50"))
 	    l.col.skat[[i]] <- rep(def.col, length(l.skat[[i]]))
	    l.col[[i]] <- rep(def.col, length(l.sing[[i]]))
            for (j in 1:length(qq.color.criteria)) {
	        crit.j <- gsub("&","&sing$",gsub("\\|","|sing$",qq.color.criteria[j]))
		l.col[[i]] <- eval(parse(text = paste0("ifelse(sing$",crit.j,", qq.color[j], l.col[[i]])")))
            }
        }
	rm(sing,skat)
}

ls()

p.sing <- unlist(l.sing)
p.skat <- unlist(l.skat)

## Added for custom colors
if (qq.color[1]=="color" | qq.color.criteria[1] != "FALSE") {
    col.sing <- unlist(l.col)
    col.skat <- unlist(l.col.skat)
}

print("Removing objects")
rm(l.sing,l.skat)
gc()

# Calculate Genomic Control lambda
# Even if using the offset, the GC lambda will be for the whole dataset
print("Calculating lambdas")
lambda.sing <- getLambda(p.sing)
lambda.skat <- getLambda(p.skat)

if (use.lambda) {
	main.sing <- paste("GC Lambda:",lambda.sing)
	main.skat <- paste("GC Lambda:",lambda.skat)
} else {
	main.sing <- main.skat <- ""
}

print("Plotting singleSNP")
# png(paste0(data.file,"_QQPLOT.png"))
eval(parse(text=paste0(output.type,"(\"qqplot_sing.",output.type,"\")")))
# png("singlesnp_qq.png")
if (qq.color[1] != "color" & qq.color.criteria[1] == "FALSE") {
    qqunif(p.sing, main=main.sing, col=qq.color, ci=include.ci)
} else {
    qqunif(p.sing, main=main.sing, col=col.sing, ci=include.ci)
}
# qqtrunc(p.sing,threshold=offset,main=main.sing,col=qq.color,ci=include.ci)
dev.off()

print("Plotting skat")
eval(parse(text=paste0(output.type,"(\"qqplot_skat.",output.type,"\")")))
# png("skat_qq.png")
if (qq.color[1] != "color" & qq.color.criteria[1] == "FALSE") {
    qqunif(p.skat, main=main.sing, col=qq.color, ci=include.ci)
} else {
    qqunif(p.skat, main=main.sing, col=col.skat, ci=include.ci)
}
# qqtrunc(p.skat,threshold=offset,main=main.skat,col=qq.color,ci=include.ci)
dev.off()





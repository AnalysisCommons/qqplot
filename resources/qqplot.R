LOCAL <- FALSE
if(LOCAL) {
  # data.files <- Sys.glob("~/DNAnexus/data/mb_test_data_chr*.Rda")
  data.files <- "~/DNAnexus/data/mb_test_data_chr22_color.txt"
  filter <- "FALSE"
  p.col <- "p"
  use.lambda <- TRUE
  include.ci <- TRUE
  qq.color <- c("red","yellow")
  qq.color.criteria <- c("maf<0.1 & maf>0.05","maf<=0.05")
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
get.delim <- function(f) {
       	d <- readLines(f,n=1)
	delim <- ifelse(grepl(" ", d)," ",ifelse(grepl(",",d),",","\t"))
        return(delim)
}

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
        
l.sing <- vector(mode="list",length=length(data.files))

## Added for custom colors
if (qq.color[1]=="color" | qq.color.criteria[1]!="FALSE") {
    l.col <- vector(mode="list",length=length(data.files))
}

# Although the data may not be singleSNP all the time, 'sing' is used to be consistent between R files
for (i in 1:length(data.files)) {
	cat("Loading ",data.files[i],"\n")
        # Read the first few lines to check the delimiter, then read the file
        delim <- get.delim(data.files[i])
        sing <- read.table(data.files[i], header=T, as.is=T, sep=delim)

        # Old code: check filetype in filename
	# if (grepl("txt$",data.files[i])) {
	# 	sing <- read.table(data.files[i],header=T,as.is=T)
	# } else {
	# 	sing <- read.csv(data.files[i],header=T,as.is=T)
	# }

	# If using a filter, subset the data here
	if (filter != "FALSE") {
          sing <- eval(parse(text = paste0("subset(sing,",filter,")")))
	}

	l.sing[[i]] <- sing[,p.col]

        ## Added for custom colors
	if (qq.color[1]=="color") {
	    l.col[[i]] <- sing[,"color"]
	} else if (qq.color.criteria[1]!="FALSE") {
	    def.col <- ifelse(!is.element("black",qq.color), "black", ifelse(!is.element("blue",qq.color),"blue","gray50"))
            l.col[[i]] <- rep(def.col, length(l.sing[[i]]))
            for (j in 1:length(qq.color.criteria)) {
		crit.j <- gsub("&","&sing$",gsub("\\|","|sing$",qq.color.criteria[j]))
                l.col[[i]] <- eval(parse(text = paste0("ifelse(sing$",crit.j,", qq.color[j], l.col[[i]])")))
            }
	}
	rm(sing)
}

ls()

p.sing <- unlist(l.sing)

## Added for custom colors
if (qq.color[1]=="color" | qq.color.criteria[1] != "FALSE") {
    col.sing <- unlist(l.col)
}


print("Removing objects")
rm(l.sing)
## Added for custom colors
if (qq.color[1]=="color") {
    rm(l.col)
}
gc()

# Calculate Genomic Control lambda
# Even if using the offset, the GC lambda will be for the whole dataset
print("Calculating lambdas")
lambda.sing <- getLambda(p.sing)

if (use.lambda) {
	main.sing <- paste("GC Lambda:",lambda.sing)
} else {
	main.sing <- ""
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





# R_SEBAL_calculate_Rn

Determine the input directory for the albedo data and output directories for Rn1030 and Rn24

```R
indir.alb = "G:/large_datasets/USA/california/modis_fluxtower_sites/MCD43A3/projected/"
indir.GMAO = "G:/large_datasets/USA/california/modis_fluxtower_sites/GMAO/"
outdir.1030 = "G:/large_datasets/USA/california/modis_fluxtower_sites/SEBAL/Rn/1030am/"
outdir.24 = "G:/large_datasets/USA/california/modis_fluxtower_sites/SEBAL/Rn/24hmean/"
```

Load the albedo list
```R
setwd(indir.alb)
file.list.alb.in = list.files(pattern="WSA")
file.list.alb = file.list.alb.in[grep("shortwave.tif",file.list.alb.in)]
yyyyjjj.albedo = as.numeric(substr(file.list.alb,10,16))
 
 # Set the start and end dates for generating Rn
start.end = c(2011061,2011180) #  Year and julian day

# Generate a subset of the filelist for just the period of interest

file.list.alb.sub = file.list.alb[(yyyyjjj.albedo>=start.end[1]) & (yyyyjjj.albedo<start.end[2])]
yyyyjjj.alb.sub = as.numeric(substr(file.list.alb.sub,10,16))

```
Load example albedo file, serves as a template for resampling MERRA data
```R
x=file.list.alb.sub[1]
balc = raster(x)
balc[balc==32767]=NA
```
Generate lists of the GMAO SW and LW files
```R
patt="nc"
flist.swlw <- list.files(indir.sw,pattern = patt)
dates.swlw = strptime(substr(flist.swlw,37,45),format="%Y%m%d")
yyyyjjj.swlw = as.numeric(format(dates.swlw,"%Y%j"))
index.swlw = which((yyyyjjj.swlw>=start.end[1])&(yyyyjjj.swlw<=start.end[2]))
file.list.swlw.sub = flist.swlw[index.swlw]  # File names in yyyymmdd
yyyyjjj.swlw.sub=yyyyjjj.swlw[index.swlw]

flist.lwup = list.files(indir.lwup,pattern="nc")
dates.lwup = strptime(substr(flist.lwup,37,44),format="%Y%m%d")
yyyyjjj.lwup = as.numeric(format(dates.lwup,"%Y%j"))
index.lwup = which((yyyyjjj.lwup>=start.end[1])&(yyyyjjj.lwup<=start.end[2]))
file.list.lwup.sub = flist.lwup[index.lwup]  # File names in yyyymmdd
yyyyjjj.lwup.sub=yyyyjjj.lwup[index.lwup]
```

Loop through all days in the list for SW radiation
```R
cnt=0
N = length(file.list.swlw.sub) # For whole list
for (x in 1:(N-1)) {
	setwd(indir.sw)
	bsw.1030am = raster(file.list.swlw.sub[x], band=19, varname = "SWGDN")  # Shortwave at ground down at 11am (PST)
	# http://globalchange.nasa.gov/KeywordSearch/Metadata.do?Portal=daacs&KeywordPath=[Keyword%3D%27LWGAB%27]&EntryId=GES_DISC_MATUNXRAD_V5.2.0&MetadataView=Text&MetadataType=0&lbnode=mdlb3

	bsw.24 = stack(file.list.swlw.sub[x],varname="SWGDN")
  bsw.24a = stack(file.list.swlw.sub[x],bands=c(8:24),varname="SWGDN")
  bsw.24b = stack(file.list.swlw.sub[x+1],bands=c(1:7),varname="SWGDN")
  bsw.ab = stack(bsw.24a,bsw.24b)
	bsw.24mean = round(mean(bsw.ab),1)  # Calculates the mean of the stack
	cnt=cnt+1
	print(paste("sw ", cnt," of ", N, sep=""))
  flush.console()  # Prints output to screen during the loop run
	projection(bsw.1030am)=projection(balc)
	projection(bsw.24mean)=projection(balc)
	bsw.1030.resam = round(resample(bsw.1030am,balc,method="bilinear"),1)
	bsw.24.resam = round(resample(bsw.24mean,balc,method="bilinear"),1)
	if (cnt==1){
		sw1030stack = stack(bsw.1030.resam)
		sw24stack = stack(bsw.24.resam)
		}
	if (cnt>1){
		sw1030stack = addLayer(sw1030stack,bsw.1030.resam)
		sw24stack = addLayer(sw24stack,bsw.24.resam)
		}
}

writeRaster(sw24stack,paste(outdir.24,"swinWm2_24_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)
writeRaster(sw1030stack,paste(outdir.1030,"swinWm2_1030_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)

rm(sw24stack)  # clears up memory
rm(sw1030stack)
```
Do the same for longwave

```R
cnt=0
for (x in 1:(N-1)) {
	# LW in
	blw.in.1030am = raster(file.list.lwup.sub[x], band=19, varname = "lwgab") # longwave absorbed = LWin
	blw.in.24 = stack(file.list.lwup.sub[x],varname="lwgab")
	blw.in.24a = stack(file.list.swlw.sub[x],bands=c(8:24),varname="LWGAB")
  blw.in.24b = stack(file.list.swlw.sub[x+1],bands=c(1:7),varname="LWGAB")
  blw.in.ab = stack(blw.in.24a,blw.in.24b)
	blw.in.24mean = round(mean(blw.in.ab),1)  # Calculates the mean of the stack
	
	# LW out
	blw.out.1030am = raster(file.list.lwup.sub[x], band=19, varname = "lwgem") #  LWGEM = longwave emitted = LWout
	blw.in.24a = stack(file.list.swlw.sub[x],bands=c(8:24),varname="LWGEM")
  blw.out.24b = stack(file.list.swlw.sub[x+1],bands=c(1:7),varname="LWGEM")
  blw.out.ab = stack(blw.in.24a,blw.in.24b)
	blw.out.24mean = round(mean(blw.out.ab),1)  # Calculates the mean of the stack
	cnt=cnt+1
	print(paste("lw ", cnt," of ", N, sep=""))
  	flush.console()
	projection(blw.in.1030am)=projection(balc)
	projection(blw.out.1030am)=projection(balc)
	projection(blw.in.24mean)=projection(balc)
	projection(blw.out.24mean)=projection(balc)
	blw.in.1030.resam = round(resample(blw.in.1030am,balc,method="bilinear"),1)
	blw.in.24.resam = round(resample(blw.in.24mean,balc,method="bilinear"),1)
	blw.out.1030.resam = round(resample(blw.out.1030am,balc,method="bilinear"),1)
	blw.out.24.resam = round(resample(blw.out.24mean,balc,method="bilinear"),1)
if (cnt==1){
	lwin1030stack = stack(blw.in.1030.resam)
	lwin24stack = stack(blw.in.24.resam)
	lwout1030stack = stack(blw.out.1030.resam)
	lwout24stack = stack(blw.out.24.resam)
	}
if (cnt>1){
	lwin1030stack = addLayer(lwin1030stack,blw.in.1030.resam)
	lwout1030stack = addLayer(lwout1030stack,blw.out.1030.resam)
	lwin24stack = addLayer(lwin24stack,blw.in.24.resam)
	lwout24stack = addLayer(lwout24stack,blw.out.24.resam)
		}
}

writeRaster(lwin24stack,paste(outdir.24,"lwinWm2_24_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)
writeRaster(lwout24stack,paste(outdir.24,"lwoutWm2_24_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)
writeRaster(lwin1030stack,paste(outdir.1030,"lwinWm2_1030_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)
writeRaster(lwout1030stack,paste(outdir.1030,"lwoutWm2_1030_",as.character(start.end[1]),"_",as.character(start.end[2]),sep=""),overwrite=TRUE)
```


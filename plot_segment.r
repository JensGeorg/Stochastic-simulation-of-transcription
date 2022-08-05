plot_segment<-function(segments){
  o<-segments
  plot(o$position,o$delay,ylim=c(0, max(o$delay)),pch=19, xlab="position [nt]", ylab="delay [min]")
  segs<-unique(o$delay_fragment)
  oo<-grep("O",segs)
  if(length(oo)>0){
    segs<-segs[-oo]
  }
  for(j in 1:length(segs)){
    temp<-which(o$delay_fragment==segs[j])
    temp<-o[temp,]
    intercept<-temp$intercept[1]
    slope<-temp$slope[1]
    s<-min(temp$position)
    e<-max(temp$position)
    curve(slope*x+intercept, add=T, from=s, to=e, lwd=2)
	mid<-s+(e-s)/2
	y_mid<-slope*mid+intercept+0.5
	velocity<-round(1/slope/60,1)
	text(mid, y_mid, labels=paste0(velocity, " nt/s"))
  }    
}

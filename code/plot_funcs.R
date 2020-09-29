export_theme <- theme_tufte() + 
  theme(
    ## Axis text and titles
    axis.text.x = element_text(size=8,family="sans"),
    axis.text.y=element_text(size=8,family="sans"),
    axis.title.x=element_text(size=8,family="sans",vjust=-1),
    axis.title.y=element_text(size=8,family="sans"),
    
    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    
    ## Title
    plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
    
    ## Legends
    legend.title=element_text(size=8,family="sans",face="italic"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans"),
    strip.background=element_rect(fill="#f0f0f0"))
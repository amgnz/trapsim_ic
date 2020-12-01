#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            TrapSim Tool - Single species - Island Conservation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Code developed by Andrew Gormley, Manaaki Whenua - Landcare Research
#Tailored for the Island conservation coati eradication project.

#Includes Audrey's grid based approach as well - Audrey says it works with cells of 500m and even 200m...(?)
#28/8 - some testing would suggest that for nightly checking, setting the grid square dimension at 4*sigma is best when compared to IBM


library("shiny")
library("shinythemes")
library("leaflet")
library("maptools")
library("spatstat") #For the owin command...and for runifpoint
library("RColorBrewer")
library("leaflet")
library("rgdal")
library("rgeos")
library("proxy")
library("sf")
library("raster")
library("DT")

# setwd("C:\\Users\\gormleya\\OneDrive - MWLR\\Documents\\CAEM\\IslandConservation\\TrapSimFeasibility\\Shiny")

getshp<-"Robinson_Coati"
shp.zones<-"Robinson_Coati_Workzones"
ras.2<-raster("habitat_specific.asc")
# ras.z2<-raster("huntingeffort.asc")


cols.vec<-brewer.pal(8,"Paired")
cols.eff<-brewer.pal(8,"Reds")
proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
#4326 is the EPSG for WGS84...

#~~~~~~~~~~~~~~~ A bunch o' functions ~~~~~~~~~~~~~~~~~~~~~~~~
#Define the resampling function
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

trap.cost.func<-function(a,b,c,d,e){
  #e.g. trap.cost.func(a=check.interval, b=n.traps, c=input$traps.per.day, d=input$day.rate,e=cost.per.trap)
  checks<-length(a)  #The number of trap checks to do = plus setup
  trap.cost<-as.integer((((b*checks)/c)*d)+(b*e))
  return(trap.cost)
}

hunt.cost.func<-function(a,b){
  day.rate<-a
  days.zone<-b
  hunt.cost<-day.rate*sum(days.zone)
  return(hunt.cost)
}

get.alpha.beta<-function(g0.mean, g0.sd){
  g0.var<-g0.sd^2
  paren<-    (g0.var + g0.mean^2 - g0.mean)
  alpha<- -g0.mean*paren/g0.var
  beta<- paren*(g0.mean-1)/g0.var
  return(list(alpha=alpha, beta=beta))
}

get.loc.shape<-function(sigma.mean, sigma.sd){
  m<-sigma.mean
  s<-sigma.sd
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  return(list(location=location, shape=shape))
}

#Function to make animal locations from a raster of relative abundance/habitat.
get.pest.locs<-function(ras.2, n.poss, shp){
  
  gridpolygon<-rasterToPolygons(ras.2)
  grids<-length(gridpolygon@data[,1])
  pest.per.grid<-round(gridpolygon@data[,1]*n.poss/grids,0)   #We want to sample more than we need becase we are going to remove some that are outside the shapefile
  coords<-vector("list", grids) #Set up to save the coordinates - cant work out how to get spsample to sample across all grids
  for (i in 1: grids){
    if(pest.per.grid[i]>0){
      coords[[i]]<-spsample(gridpolygon@polygons[[i]], n =pest.per.grid[i] , type = "random")@coords
    }
  }
  coords<-as.data.frame(do.call("rbind", coords)) #Put them altogether from the list
  coords<-coords[inside.owin(coords[,1], coords[,2], shp),]  #Remove the ones from outside the shapefile...
  coords<-coords[sample(1:dim(coords)[1], size=n.poss, replace=FALSE),]  #Then sample to get the desired actual sample size
  
  return(coords)
  
}


#  Make the trap locations from a spacing, buffer and shapfile
make.trap.locs<-function(x.space,y.space,buff,shp){
  
  b.box<-bbox(shp)
  
  traps.x<-seq(from=b.box[1,1], by=x.space, to=b.box[1,2])
  traps.y<-seq(from=b.box[2,1], by=y.space, to=b.box[2,2])
  
  traps<-as.data.frame((expand.grid(traps.x, traps.y)))
  colnames(traps)<-c("X","Y")
  shp.buff<-gBuffer(shp,width=-buff)
  #Remove traps that are outside the window,,,
  traps<-traps[inside.owin(traps[,1], traps[,2], shp.buff),]
  
  return(traps)
}


# Function to split the input boxes for scenarios
split.str<-function(a){
  b<-as.numeric((strsplit(as.character(a),split= "/"))[[1]])
  return(b)
}
#~~~~~~~~~~~~~~~~~  End functions ~~~~~~~~~~~~~~~~~~~~~~


server<-function(input, output, session) {
  
  # mydata.ini<-reactive({
  #   #~~~~~~~~~~~~~DEFINE THE AREA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   #~~~~~~~~ Working with a area of specifed size, or a shapefile...
  #   # validate(
  #   #   need(input$area.ha != "", "Please enter a size for the area in hectares")
  #   # )
  #   
  #   
  #   # if(input$area_type=="Area"){
  #   #   ha.tmp<-input$area.ha
  #   #   size<-sqrt(ha.tmp*10000)
  #   #   xy<-c(1730000,5620000)
  #   #   a<-cbind(c(0,size,size,0, 0)+xy[1],c(0,0,size,size,0)+xy[2])
  #   #   ID<-"a"
  #   #   shp<-SpatialPolygons(list(Polygons(list(Polygon(a)),ID)), proj4string=CRS(proj4string))
  #   #   
  #   # }else{
  #     #Read in the map - but start with a tmp size to get things happening...
  #     
  #     # if(is.null(input$shp.file)){
  #     #   ha.tmp<-input$area.ha
  #     #   size<-sqrt(ha.tmp*10000)
  #     #   xy<-c(1730000,5620000)
  #     #   a<-cbind(c(0,size,size,0, 0)+xy[1],c(0,0,size,size,0)+xy[2])
  #     #   ID<-"a"
  #     #   shp<-SpatialPolygons(list(Polygons(list(Polygon(a)),ID)), proj4string=CRS(proj4string))
  #     #   
  #     # }else{
  #     #   myshape<-input$shp.file
  #     #   dir<-dirname(myshape[1,4])
  #     #   for ( i in 1:nrow(myshape)) {
  #     #     file.rename(myshape[i,4], paste0(dir,"/",myshape[i,1]))
  #     #   }
  #     #   getshp <- list.files(dir, pattern="*.shp", full.names=TRUE)
  #     #   shp<-readOGR(getshp)
  #     #   proj4string<-crs(shp)  #get the proj string from the shapefile itself
  #     #   shp<-gBuffer(shp, width=1) #Can fix some orphaned holes issues
  #     #   # shp<-shp[2,]
  #     # }
  #     #Temporary whilst above is commented out      
  #     shp<-readOGR("Shapefiles",getshp)
  #     proj4string<-crs(shp)  #get the proj string from the shapefile itself
  #     shp<-gBuffer(shp, width=1) #Can fix some orphaned holes issues
  #     
  #     # shp.2<-readOGR("Shapefiles",shp.zones)
  #     # proj4string<-crs(shp.2)  #get the proj string from the shapefile itself
  #     # shp.2<-gBuffer(shp.2, width=1) #Can fix some orphaned holes issues
  #     
  #   # }
  #   
  #   
  #   b.box<-bbox(shp)
  #   ha<-sapply(slot(shp, "polygons"), slot, "area")/10000  
  #   
  #   
  #   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   
  #   traps.x.space.vec<-split.str(input$traps.x.space)
  #   traps.x.space<-traps.x.space.vec[1]
  #   # traps.x.space<-as.numeric((strsplit(as.character(traps.x.space.tmp),split= "/"))[[1]])[1]
  #   traps.y.space.vec<-split.str(input$traps.y.space)
  #   traps.y.space<-traps.y.space.vec[1]
  #   # traps.y.space<-input$traps.y.space
  #   buff<-input$traps.buff
  #   
  #   
  #   validate(
  #     need(input$traps.buff != "", "Please enter a value for the trap buffer"),
  #     need(input$traps.x.space != "", "Please enter a value for the X trap spacing"),
  #     need(input$traps.y.space != "", "Please enter a value for the Y trap spacing"),
  #     need(input$n.nights != "", "Please enter a value for the number of Nights"),
  #     need(input$n.check != "", "Please enter a value for the Check interval"),
  #     need(input$p.bycatch != "", "Please enter a value for the Daily bycatch rate"),
  #     need(input$max.catch != "", "Please enter a value for the prob of Max catch")
  #   )
  #   
  #   if(input$trap_type=="Sim"){
  #     #~~~~Now make the traps from the specified spacings...
  #     #This is for plotting on the map.i.e. all traps and grid cell...
  #     #Make the traps across the entire bounding box
  #     traps.x<-seq(from=b.box[1,1], by=traps.x.space, to=b.box[1,2])
  #     traps.y<-seq(from=b.box[2,1], by=traps.y.space, to=b.box[2,2])
  #     
  #     traps<-as.data.frame((expand.grid(traps.x, traps.y)))
  #     colnames(traps)<-c("X","Y")
  #     
  #   }else{
  #     
  #     if(is.null(input$trap.file)){
  #       x<-c(1730000)
  #       y<-c(5620000)
  #       
  #       traps<-data.frame("X"=x, "Y"=y)
  #       colnames(traps)<-c("X","Y")
  #       
  #     }else{
  #       infile <- input$trap.file
  #       data.in<-read.csv(infile$datapath) 
  #       idx<-match(c("X", "Y"),toupper(colnames(data.in)))
  #       if(is.na(idx[1])==TRUE){
  #         idx<-match(c("EASTING","NORTHING"),toupper(colnames(data.in)))
  #       }
  #       traps<-data.in[,idx]
  #       colnames(traps)<-c("X","Y")
  #     }
  #   }
  #   
  #   
  #   shp.buff<-gBuffer(shp,width=-buff)
  #   #Remove traps that are outside the window,,,
  #   traps<-traps[inside.owin(traps[,1], traps[,2], shp.buff),]
  #   
  #   #Remove traps that are outside the window - need this for irregular shapefiles.
  #   # traps<-traps[inside.owin(traps[,1], traps[,2], shp),]
  #   #How many traps are left.
  #   n.traps<-dim(traps)[1]   
  #   
  #   
  #   
  #   #Temporary pest animal coordinates...
  #   n.poss<-input$numb.poss
  #   
  #   
  #   #Temporary commented out.    
  #   # if(is.null(input$ras.1)==TRUE){
  #   if((input$ras.1)=="Random"){
  #     #1. Random locations.
  #     n.poss.tmp<-(runifpoint(n.poss,shp))
  #     animals.xy.ini<-as.data.frame(n.poss.tmp)
  #   }else{
  #     #2. Grid specific densities
  #     # infile.ras <-input$ras.1
  #     # ras.2<-raster(infile.ras$datapath)
  #     #Call the function
  #     animals.xy.ini<-get.pest.locs(ras.2, n.poss, shp)
  #   }
  #   
  #   colnames(animals.xy.ini)<-c("X","Y")
  #   
  #   
  #   
  #   
  #   # #How long to run the simulation for.
  #   # trap.nights<-split.str(input$trap.nights)[1]
  #   # trap.start<-split.str(input$trap.start)[1]
  #   # #When are traps checked 
  #   # n.check<-split.str(input$n.check)[1]
  #   # #This sets the trap checking interval. i.e. traps are cleared and reset on these nights only...
  #   # check.interval<-seq(from=trap.start, to=(trap.start+trap.nights), by=n.check)
  #   # 
  #   # trap.cost<-trap.cost.func(a=check.interval,b=n.traps,c=input$traps.per.day, d=input$day.rate, e=input$cost.per.trap)
  #   # #e.g. trap.cost.func(a=check.interval, b=n.traps, c=input$traps.per.day, d=input$day.rate,)
  #   
  # 
  #   return(list(traps=traps, shp=shp, ha=ha, animals.xy.ini=animals.xy.ini, p4s=proj4string))      
  # })
  
  
  #Read in the shapefile  
  mydata.shp<-reactive({
    shp<-readOGR("Shapefiles",getshp)
    proj4string<-crs(shp)  #get the proj string from the shapefile itself
    shp<-gBuffer(shp, width=1) #Can fix some orphaned holes issues
    ha<-sapply(slot(shp, "polygons"), slot, "area")/10000  
    return(list(shp=shp, p4s=proj4string, ha=ha))
  })
  
  
  #This sets up the zones for hunting effort.  
  mydata.zone<-reactive({
    shp.2<-readOGR("Shapefiles",shp.zones)
    proj4string<-crs(shp.2)  #get the proj string from the shapefile itself
    
    effort.zone<-c(input$effort.a, input$effort.b, input$effort.c, input$effort.d)
    theta.hat<-1-exp(-((input$hunt.rho*log(effort.zone))^input$hunt.k))
    p.hunt.day<-theta.hat
    
    return(list(shp.2=shp.2, p.hunt.day=p.hunt.day))    
  })
  
  #Make the table of scenarios 
  mydata.scen<-reactive({
    #We nmeed to call in the inputs and then make the params... and run it over all of those.
    #~~~~~~~~~~~~~~~~~~~~~~~~TRAPS~~~~~~~~~~~~~~~~~~
    #Trap spacing...
    traps.x.space.vec<-as.integer(split.str(input$traps.x.space))  #split str uses the "/" for splitting
    traps.y.space.vec<-as.integer(split.str(input$traps.y.space))
    #Trapping period
    trap.start.vec<-as.integer(split.str(input$trap.start))
    trap.nights.vec<-as.integer(split.str(input$trap.nights))
    n.check.vec<-as.integer(split.str(input$n.check))
    
    hunt.days.a.vec<-as.integer(split.str(input$hunt.days.a))
    
    params<-expand.grid(traps.x.space.vec,traps.y.space.vec, trap.start.vec, trap.nights.vec, n.check.vec, hunt.days.a.vec)
    names(params)<-c("x.space", "y.space", "trap.start","trap.nights","check.interval", "hunt.days.a")
    
    params$Scenario<-rownames(params)
    params<-params[,c(ncol(params),1:(ncol(params)-1))]
    
    return(list(params=params))
  })
  
  
  # Make the traps and animals on the interactive map. Not linked to the actual simulation
  mydata.map<-reactive({
    shp<-mydata.shp()$shp    
    
    traps.x.space<-as.numeric(input$traps.x.space.i)
    traps.y.space<-as.numeric(input$traps.y.space.i)
    # traps.y.space<-input$traps.y.space
    buff<-as.numeric(input$traps.buff.i)
    traps<-make.trap.locs(traps.x.space, traps.y.space, buff, shp)
    
    #Temporary pest animal coordinates...
    n.poss<-input$numb.poss.i
    #Temporary commented out.    
    # if(is.null(input$ras.1)==TRUE){
    if((input$ras.1)=="Random"){
      #1. Random locations.
      n.poss.tmp<-(runifpoint(n.poss,shp))
      animals.xy.ini<-as.data.frame(n.poss.tmp)
    }else{
      #2. Grid specific densities
      # infile.ras <-input$ras.1
      # ras.2<-raster(infile.ras$datapath)
      #Call the function
      animals.xy.ini<-get.pest.locs(ras.2, n.poss, shp)
    }
    colnames(animals.xy.ini)<-c("X","Y")
    
    return(list(traps=traps, animals.xy.ini=animals.xy.ini))
  })
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~  Now simulate the actual trapping...~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  datab<-eventReactive(input$act.btn.trapsim,{
    
 #Some validation stuff
    validate(
      need(input$traps.buff != "", "Please enter a value for the trap buffer"),
      need(input$traps.x.space != "", "Please enter a value for the X trap spacing"),
      need(input$traps.y.space != "", "Please enter a value for the Y trap spacing"),
      need(input$n.nights != "", "Please enter a value for the number of Nights"),
      need(input$n.check != "", "Please enter a value for the Check interval"),
      need(input$p.bycatch != "", "Please enter a value for the Daily bycatch rate"),
      need(input$max.catch != "", "Please enter a value for the prob of Max catch")
    )
    
    #Get the parameter values for the scenarios
    params<-mydata.scen()$params
    n.scen<-dim(params)[1]    
    pop.size.list<-vector("list",n.scen)
    pop.zone.list<-vector("list",n.scen)
    trap.catch.list<-vector("list",n.scen)
    hunt.catch.list<-vector("list",n.scen)
    # pop.size.zone.mat[[ii]][,t+1]
    
    withProgress(message="Running simulation ",value=0,{
      
      for(kk in 1:n.scen){
        incProgress(kk/n.scen, detail = paste("Doing scenario ", kk," of", n.scen))
        
        # #Trapping period
        trap.start<-params$trap.start[kk]
        trap.nights<-params$trap.nights[kk]
        
        trap.period<-seq(from=trap.start, to=(trap.start+trap.nights-1), by=1)
        
        # #When are traps checked 
        n.check<-params$check.interval[kk]
        
        
        #How long to run the simulation for.
        n.nights<-input$n.nights
        
        # #This sets the trap checking interval. i.e. traps are cleared and reset on these nights only...
        check.interval<-seq(from=trap.start, to=(trap.start+trap.nights), by=n.check)
        p.bycatch<-input$p.bycatch
        
        
        #Probability of being hunted per day
        p.hunt.day<-mydata.zone()$p.hunt.day
        
        #What days are the animals being hunted??
        days.zone<-c(params$hunt.days.a[kk], input$hunt.days.b, input$hunt.days.c, input$hunt.days.d)
        
        
        hunt.cost.sim<-hunt.cost.func(input$day.rate.hunt, days.zone)
        
        
        hunt.periods<-vector("list",4)
        hunt.periods[[1]]<-seq(from=input$hunt.start.a, (input$hunt.start.a+days.zone[1]-1))
        hunt.periods[[2]]<-seq(from=input$hunt.start.b, (input$hunt.start.b+days.zone[2]-1))
        hunt.periods[[3]]<-seq(from=input$hunt.start.c, (input$hunt.start.c+days.zone[3]-1))
        hunt.periods[[4]]<-seq(from=input$hunt.start.d, (input$hunt.start.d+days.zone[4]-1))
        
        
        
        g0.mean<-input$g0.mean
        g0.sd<-input$g0.sd
        sigma.mean<-input$sigma.mean 
        sigma.sd<-input$sigma.sd
        rmax.poss<-input$rmax.poss
        K.poss<-10 #input$K.poss
        #When the reproductive period starts and how long it lasts (i.e. spread out over...)
        rep.start<-input$rep.start
        rep.nights<-input$rep.nights
        rep.tmp<-seq(from=rep.start, by=1, to=(rep.nights+rep.start-1))
        aa<-ceiling(n.nights/365) #basically the number of years
        bb<-((1:aa)-1)*365
        repro.interval<-rep(rep.tmp,aa)+rep(bb, each=rep.nights)
        rep.start.vec<-(bb+rep.start)
        
        ha<-mydata.shp()$ha
        shp<-mydata.shp()$shp
        shp.2<-mydata.zone()$shp.2
        
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make the trap locations
        traps<-make.trap.locs(params$x.space[kk], params$y.space[kk],input$traps.buff,shp)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        #Carrying capacity for the area
        K.tot<-K.poss*ha
        
        # 
        # if(input$animal_type=='Dens'){
        #   dens.poss<-input$dens.poss
        #   #~~~ Make the animals at the required density.
        #   n.poss<-floor(dens.poss*ha)
        #   
        # }else{
        #   n.poss<-input$numb.poss
        # }
        
        n.poss<-input$numb.poss
        
        max.catch<-input$max.catch
        # traps<-mydata()$traps
        coordinates(traps) <- c( "X", "Y" )
        proj4string(traps) <- CRS(proj4string)
        traps.xy<-as.data.frame(traps)
        n.traps<-dim(traps.xy)[1]
        
        
        trap.cost.sim<-trap.cost.func(a=check.interval,b=n.traps,c=input$traps.per.day, d=input$day.rate, e=input$cost.per.trap)
        
        
        
        if (input$sim_type=='grid'){
          cell.width<-input$cell.width
          cell.area.m2<-cell.width^2
          
          #Make the grid...
          shp<-st_as_sf(shp)
          shp_grid<-shp%>%st_make_grid(cellsize=cell.width, what="polygons")%>%st_intersection(shp)
          grid.traps.master<-lengths(st_intersects(shp_grid, st_as_sf(traps)))#*max.catch  #Number in each cell
          
        }
        
        n_its<-input$n.its
        pop.size.mat<-matrix(NA,nrow=n_its,ncol=n.nights+1)
        trap.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of trapping - for each iteration - how many that night/day
        hunt.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of hunting
        pop.size.zone.vec<-vector("list",n_its)
        
        
        for(ii in 1:n_its){
          pop.size.mat[ii,1]<-n.poss
          
          #~~~~~~~~~Make some animals~~~~~~~~~~
          # if(is.null(input$ras.1)==TRUE){
          if((input$ras.1)=="Random"){
            #1. Random locations.
            n.poss.tmp<-(runifpoint(n.poss,shp))
            animals.xy<-as.data.frame(n.poss.tmp)
          }else{
            #2. Grid specific densities
            # infile.ras <-input$ras.1
            # ras.2<-raster(infile.ras$datapath)
            #Call the function
            animals.xy<-get.pest.locs(ras.2, n.poss, shp)
          }
          
          
          colnames(animals.xy)<-c("X","Y")
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          # n.poss.tmp<-(runifpoint(n.poss,shp))
          # animals.xy<-as.data.frame(n.poss.tmp)
          # colnames(animals.xy)<-c("X","Y")
          n.animals<-dim(animals.xy)[1]
          animals.xy$Dead<-0
          coordinates(animals.xy) <- ~ X+Y
          proj4string(animals.xy) <- proj4string(shp)
          
          
          
          
          
          alpbet<-get.alpha.beta(g0.mean, g0.sd)
          animals.xy$g0u<-rbeta(n.animals, alpbet$alpha, alpbet$beta)
          animals.xy$g0u[animals.xy$g0u<0]<-0 #Probably not needed
          
          locshp<-get.loc.shape(sigma.mean, sigma.sd)
          animals.xy$Sigma<-rlnorm(n.animals ,meanlog=locshp$location, sdlog=locshp$shape)
          
          #The first one initialises for a grid based simulation - i.e. Audrey's model.
          #The seond is the trap+animal pairwise 
          if (input$sim_type=='grid'){
            animals.SP<-animals.xy
            coordinates(animals.SP) <- c( "X", "Y" )
            proj4string(animals.SP) <- CRS(proj4string)
            which.grid<-st_within(st_as_sf(animals.SP), shp_grid)
            animals.xy$CellIndex<-as.data.frame(which.grid)$col.id
            animals.xy$PreProb<-(2*pi*animals.xy$g0u*animals.xy$Sigma^2)/cell.area.m2   #The pre probability...
            grid.traps<-grid.traps.master
            grid.traps.capacity<-grid.traps.master*max.catch
          }else{
            
            dist2.xy<-matrix(NA,n.traps,n.animals)
            prob.xy<-matrix(0,n.traps,n.animals)
            dist.xy<-dist(as.data.frame(traps.xy), as.data.frame(animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
            prob.xy<-exp(-(dist.xy^2)/(2*animals.xy$Sigma^2))*animals.xy$g0u #Use the g0u for sampling...
            rm(dist.xy)
          }
          
          
          #~~~Now run the simulation of trapping the animals...
          # withProgress(message="Running Simulation",value=0,{
          trap.catch.vec<-rep(0, n.nights)
          trap.catch<-matrix(0,n.traps,n.nights)  			#Trap.catch stores the captures in each trap each night
          # trap.remain<-rep(T,n.traps)						        #Record if the trap remains in operation...TRUE, FALSE
          trap.remain<-rep(max.catch,n.traps)
          
          
          pop.size.zone.vec[[ii]]<-matrix(NA,nrow=4,ncol=n.nights+1)
          pop.size.zone.vec[[ii]][,1]<-table(over(animals.xy[animals.xy$Dead==0,], shp.2))
          
          
          for (t in 1:n.nights){								#For each night
            not.caught<-(1:n.animals)[animals.xy$Dead==0]			#Animals not already caught
            
            
            if(t%in%trap.period==TRUE){
              if(t%in%check.interval==TRUE){#If it is a trap clearance day...then reset the traps to T *before* trappig starts!
                # trap.remain<-rep(T,n.traps)		
                trap.remain<-rep(max.catch,n.traps)
                if (input$sim_type=='grid'){
                  grid.traps<-grid.traps.master
                }
                
              }
              #Turn off some of the traps according to the random probability
              # trap.remain[rbinom(n.traps,1, p.bycatch)==1]<-FALSE #Nedd to only turn off those that are on...Might be okay...
              trap.remain<-trap.remain-rbinom(n.traps,trap.remain, p.bycatch) #This modifcation deals with multiple capture traps 
              trap.remain[trap.remain<0]<-0
              
              if(sum(not.caught)>0){
                for (j in not.caught){ 							#For each animal not already caught
                  if (input$sim_type=='grid'){
                    #Based on Audrey's model...
                    pcap<-1-exp(-(animals.xy$PreProb[j]*grid.traps[animals.xy$CellIndex[j]]))
                    if(rbinom(1,1,prob=pcap)==1){ #If animal gets caught
                      animals.xy$Dead[j]<-1
                      grid.traps[animals.xy$CellIndex[j]]<-grid.traps[animals.xy$CellIndex[j]]-1
                      trap.catch.vec[t]<-trap.catch.vec[t]+1
                    }
                  }else{
                    #New code based on Multinomial (from Dean...
                    prob.tmp<-prob.xy[,j]*(trap.remain>0)				#Adjust the probabilities with trap.remain so previous traps cannot catch anything
                    cumulative.capture.prob<-1-prod(1-prob.tmp) #Cumulative probability of possum i getting captured
                    
                    if(rbinom(1,1,prob=cumulative.capture.prob)==1){#If the animal is going to get caught...
                      trap.id<-match(1,rmultinom(1,1,prob.tmp))#Which was the successful trap...
                      trap.catch[trap.id,t]<-trap.catch[trap.id,t]+1
                      # trap.animals[j]<-1					#Set trapped animal to 1
                      animals.xy$Dead[j]<-1
                      # animals.xy$Day2[j]<-t
                      # trap.remain[trap.id]<-(trap.catch[trap.id,t]<max.catch)			#Calculate whether the trap is full or not!!
                      trap.remain[trap.id]<-trap.remain[trap.id]-1
                    }
                  }
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Module for immigration to go in here...
              }
            } #End of if t %in% trap.period
            
            #Hunting
            
            
            
            hunt.vec<-c(0,0,0,0)
            zones<-c("A","B","C","D")
            for (z in 1:4){
              if(days.zone[z]>0){
                if(t%in% hunt.periods[[z]] ){
                  # idx.animal.zone<-inside.owin(animals.xy[,1], animals.xy[,2], shp.2[shp.2$WorkZone=='C',])
                  idx.animal.zone<-(inside.owin(as.data.frame(animals.xy)[,1], as.data.frame(animals.xy)[,2], shp.2[shp.2$WorkZone==zones[z],])&animals.xy$Dead==0)
                  N.tmp<-sum(idx.animal.zone)
                  if(N.tmp>0){
                    n.kill<-rbinom(n=1, p=p.hunt.day[z], size=N.tmp)
                    hunt.vec[z]<-n.kill
                    if(n.kill>0){
                      idx.dead<-sample(which(idx.animal.zone,TRUE), size=n.kill, replace=FALSE) 
                      animals.xy$Dead[idx.dead]<-1
                    }
                  }
                }
              }
            }
            hunt.catch.mat[ii,t]<-sum(hunt.vec)
            
            
            #~~~~~~~~~~~~~~~~~~~~~Reproduction~~~~~~~~~~~~~~~~~~~~~~~~~
            #Get the repdirctive vector when we hit the start of the interval...
            #If t is one of the start days of the breeding season...
            if(t%in%rep.start.vec){
              #Discrete version of rmax
              rd<-exp(rmax.poss)-1
              N0<-sum(animals.xy$Dead==0) #Current population
              K<-K.tot
              new.animals<-rd*N0*((K-N0)/K)#Number of new animals
              mu.animals<-new.animals/rep.nights #new animals per night
              if(mu.animals>0){ #Draw animals when mu is positive
                new.animal.vec<-(rpois(rep.nights, lambda=mu.animals))
                new.animal.vec[(N0+cumsum(new.animal.vec))>(K*1.1)]<-0
              }else{
                new.animal.vec<-rep(0,rep.nights) #new animals is 0 when mu is negative
              }
              yr<-match(t, rep.start.vec) #What year are we in...?
            }
            
            
            if(t%in%repro.interval==TRUE){
              #Get the index from the number of new animals
              t.idx<-t-(365*(yr-1)) #TO adjust for multi years
              N.new<-new.animal.vec[match(t.idx,repro.interval)]
              if(N.new>0){  
                #Make the new locations randomly...
                
                # if(is.null(input$ras.1)==TRUE){
                if((input$ras.1)=="Random"){
                  
                  #1. Random locations.
                  
                  new.animals.xy<-as.data.frame(runifpoint(N.new,shp))
                }else{
                  #2. Grid specific densities
                  # infile.ras <-input$ras.1
                  # ras.2<-raster(infile.ras$datapath)
                  #Call the function
                  new.animals.xy<-get.pest.locs(ras.2, N.new, shp)
                }
                
                # new.animals.xy<-as.data.frame(runifpoint(N.new,shp))
                colnames(new.animals.xy)<-c("X","Y")
                new.animals.SP<-new.animals.xy  #Why dis?
                coordinates(new.animals.SP) <- c( "X", "Y" )
                proj4string(new.animals.SP) <- CRS(proj4string)
                # new.animals.xy$g0<-g0.mean
                # new.animals.xy$Sigma<-sigma.mean
                new.animals.xy$Dead<-0
                new.animals.xy$g0u<-rbeta(N.new, alpbet$alpha, alpbet$beta)
                new.animals.xy$g0u[new.animals.xy$g0u<0]<-0.000001
                new.animals.xy$Sigma<-rlnorm(N.new ,meanlog=locshp$location, sdlog=locshp$shape)        
                
                if (input$sim_type=='grid'){
                  new.which.grid<-st_within(st_as_sf(new.animals.SP), shp_grid)
                  new.animals.xy$CellIndex<-as.data.frame(new.which.grid)$col.id
                  new.animals.xy$PreProb<-(2*pi*new.animals.xy$g0u*new.animals.xy$Sigma^2)/cell.area.m2   #The pre probability...
                }else{
                  # #Calculate the probabilities for these new ones...
                  dist.xy.new<-matrix(NA,n.traps,N.new)
                  prob.xy.new<-matrix(0,n.traps,N.new)
                  dist.xy.new<-dist(as.data.frame(traps.xy), as.data.frame(new.animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
                  prob.xy.new<-exp(-(dist.xy.new^2)/(2*new.animals.xy$Sigma^2))*new.animals.xy$g0u #Use the g0u for sampling...
                  rm(dist.xy.new)
                  prob.xy<-cbind(prob.xy, prob.xy.new)
                }

                coordinates(new.animals.xy) <- ~ X+Y
                proj4string(new.animals.xy) <- proj4string(shp)
                animals.xy<-rbind(animals.xy,new.animals.xy)            
                
                n.animals<-dim(animals.xy)[1]
              }}
            
            #Update the age-class...
            pop.size.mat[ii,t+1]<-sum(animals.xy$Dead==0)
            #This calculates the number of alive animals in each zone.  
            pop.size.zone.vec[[ii]][,t+1]<-table(over(animals.xy[animals.xy$Dead==0,], shp.2))
            
            
            
          }	 #End of the night
          if (input$sim_type=='grid'){
            trap.catch.mat[ii,]<-(trap.catch.vec)
          }else{
            trap.catch.mat[ii,]<-colSums(trap.catch)
          }
        }#End of iteration ii
        
        
        pop.zone.list[[kk]]<-apply(simplify2array(pop.size.zone.vec),c(1,2), mean)
        params$TrapCost[kk]<-trap.cost.sim
        params$HuntCost[kk]<-hunt.cost.sim
        params$TotalCost[kk]<-hunt.cost.sim + trap.cost.sim
        
        params$MeanPopSize[kk]<-round(mean(pop.size.mat[,n.nights+1]),2)
        
        pop.size.list[[kk]]<-pop.size.mat
        hunt.catch.list[[kk]]<-hunt.catch.mat
        trap.catch.list[[kk]]<-trap.catch.mat
      } #End kk    
      
      
      
    })      
    
    return(list(trap.catch=trap.catch, trap.catch.mat=trap.catch.mat, pop.size.mat=pop.size.mat, animals.xy=animals.xy, hunt.catch.mat=hunt.catch.mat, params=params, pop.size.list=pop.size.list, trap.catch.list=trap.catch.list, hunt.catch.list=hunt.catch.list, pop.zone.list=pop.zone.list))#, animals.done.xy=animals.xy))    
  })
  
  
  #~~~~~~~~~~~~~~Draw the map and plots etc ~~~~~~~~~~~~~~~~~~~~~~
  # output$plot1<-renderPlot({
  #   shp<-mydata()$shp
  #   ha<-mydata()$ha
  #   traps<-mydata()$traps
  #   animals.xy<-mydata()$animals.xy
  #   
  #   
  #   par(mar=c(2,2,7,1))
  #   plot(shp, col=cols.vec[3])
  #   if(input$area_type=="Area"){
  #     mtext("Indicative Area" ,3, cex=1.6, line=5.5) #Shapefile name
  #   }else{
  #     mtext(shp[[1]] ,3, cex=1.6, line=5.5) #Shapefile name
  #   }
  #   mtext(paste("Area =",round(ha,0), " ha"),1, line=1, cex=1.2) #Density
  #   
  #   
  #   
  #   
  #   points(traps, pch="+", cex=1, col=cols.vec[6])      
  #   points(animals.xy$X, animals.xy$Y, pch=16, col=cols.vec[2])
  #   mtext(paste("Number of traps = ",dim(traps)[1],sep=""),3, line=4, adj=0, cex=1.2)
  #   mtext(paste("Trap density = ", round(dim(traps)[1]/ha,2),"/ha",sep=""),3, line=4, adj=1, cex=1.2)
  # })
  
  
  
  
  output$mymap<-renderLeaflet({
    # shp<-mydata.shp()$shp
    shp<-mydata.zone()$shp.2
    
    shp.proj<-spTransform(shp,CRS("+proj=longlat +datum=WGS84"))
    # ha<-mydata()$ha
    
    m<-leaflet() %>%
      addTiles(group="Default")%>%
      addProviderTiles("Esri.WorldTopoMap", group = "Topo")%>%
      addProviderTiles("Esri.WorldImagery", group = "Aerial")%>%
      addPolygons(data=shp.proj, weight=2, fillColor="grey40", color="black", fillOpacity=0.2)
    m
    
  })
  
  
  
  observe({
    traps<-mydata.map()$traps
    animals.xy<-mydata.map()$animals.xy.ini
    proj4string<-mydata.shp()$p4s
    
    traps.proj<-proj4::project(traps, proj=proj4string, inverse=T) 
    tmp<-(proj4::project(animals.xy[,1:2], proj=proj4string, inverse=T))
    animals.xy$Lat<-tmp$y
    animals.xy$Lon<-tmp$x
    map<-leafletProxy("mymap")
    map%>%clearMarkers()
    
    map%>%addCircleMarkers(lng=traps.proj$x,lat=traps.proj$y, radius=3, color="black", weight=1, fill=TRUE, fillColor="red", fillOpacity=1, stroke=TRUE, group="Traps")
    map%>%addCircleMarkers(lng=animals.xy$Lon,lat=animals.xy$Lat, radius=5, color="black", weight=1, fill=TRUE, fillColor=cols.vec[2], fillOpacity=1, stroke=TRUE, group="Animals")
    
    map%>%addLayersControl(
      baseGroups = c("Default","Topo","Aerial"),
      overlayGroups=c("Traps","Animals"),
      options = layersControlOptions(collapsed = FALSE)
    )
    
  })
  
  
  output$plot2<-renderPlot({
    trap.catch.list<-datab()$trap.catch.list
    trap.catch.mat<-trap.catch.list[[as.numeric(input$result_scenario)]]
    
    # trap.catch.mat<-datab()$trap.catch.mat
    nights.vec<-1:input$n.nights
    ymax<-max(cumsum(colMeans(trap.catch.mat)))
    
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,ymax), type='n', xlab="Nights", ylab="Cumulative captures", las=1)
    for(i in 1:input$n.its){
      lines(nights.vec,cumsum(trap.catch.mat[i,]), col="grey")
    }
    lines(nights.vec,cumsum(colMeans(trap.catch.mat)), col="black")
    points(nights.vec,cumsum(colMeans(trap.catch.mat)), bg="black", pch=21)
    mtext("Cumulative captures: Trapping",3, cex=1.5, line=1)
    # mtext("Trapped per night",3, cex=1.5, line=1)
  })
  

  output$plot2.hunt<-renderPlot({
    
    hunt.catch.list<-datab()$hunt.catch.list
    trap.catch.mat<-hunt.catch.list[[as.numeric(input$result_scenario)]]
    # trap.catch.mat<-datab()$hunt.catch.mat
    nights.vec<-1:input$n.nights
    ymax<-max(cumsum(colMeans(trap.catch.mat)))
    
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,ymax), type='n', xlab="Nights", ylab="Cumulative captures", las=1)
    for(i in 1:input$n.its){
      lines(nights.vec,cumsum(trap.catch.mat[i,]), col="grey")
    }
    lines(nights.vec,cumsum(colMeans(trap.catch.mat)), col="black")
    points(nights.vec,cumsum(colMeans(trap.catch.mat)), bg="black", pch=21)
    mtext("Cumulative captures: Hunting",3, cex=1.5, line=1)
    # mtext("Trapped per night",3, cex=1.5, line=1)
  })
  
  
  
  # output$plot3<-renderPlot({
  #   trap.catch<-datab()$trap.catch
  #   par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
  #   x<-rowSums(trap.catch)
  #   max.x<-max(x)
  #   barplot(table(factor(x, levels=0:max.x)), las=1, col=cols.vec[7])
  #   mtext("Captures per trap",3, cex=1.5, line=1)
  #   
  # })
  
  observe({
    updateSelectInput(session=session, inputId="result_scenario",choices=mydata.scen()$params$Scenario)
    
  })
  
  
  output$plot4<-renderPlot({
    
    pop.size.list<-datab()$pop.size.list
    pop.size.mat<-pop.size.list[[as.numeric(input$result_scenario)]]
    # pop.size.mat<-datab()$pop.size.mat
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,max(pop.size.mat)), type='n', xlab="Nights", ylab="Population Size", las=1)
    grid()
    for(i in 1:input$n.its){
      lines(0:input$n.nights,pop.size.mat[i,], col="grey")
    }
    lines(0:input$n.nights,colMeans(pop.size.mat), col="black")
    points(0:input$n.nights,colMeans(pop.size.mat), bg="black", pch=21)
    mtext("Population size",3, cex=1.5, line=1)
    
  })
  
  
  
  
  output$plotzone<-renderPlot({
    pop.zone.list<-datab()$pop.zone.list
    pop.size.mat<-pop.zone.list[[as.numeric(input$result_scenario)]]
    # pop.size.mat<-datab()$pop.size.mat
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,max(pop.size.mat)), type='n', xlab="Nights", ylab="Population Size by Zone", las=1)
    grid()
    lines(0:input$n.nights,pop.size.mat[1,], col="red")
    lines(0:input$n.nights,pop.size.mat[2,], col="green")
    lines(0:input$n.nights,pop.size.mat[3,], col="blue")
    lines(0:input$n.nights,pop.size.mat[4,], col="purple")
    mtext("Population size",3, cex=1.5, line=1)
    legend("topright", col=c("red","green","blue","purple"), legend=c("A","B","C","D"), bty="n", lty=1)
    
  }) 
  
  # # Plot of captures per night...
  # output$plot5<-renderPlot({
  #   trap.catch.mat<-datab()$trap.catch.mat
  #   # traps<-mydata()$traps
  #   # n.traps<-dim(traps)[1]
  #   nights.vec<-1:input$n.nights
  #   ymax<-max(colMeans(trap.catch.mat))
  #   
  #   par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
  #   plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,ymax), type='n', xlab="Nights", ylab="Trapped", las=1)
  #   for(i in 1:input$n.its){
  #     lines(nights.vec,trap.catch.mat[i,], col="grey")
  #   }
  #   lines(nights.vec,colMeans(trap.catch.mat), col="black")
  #   points(nights.vec,colMeans(trap.catch.mat), bg="black", pch=21)
  #   # mtext("Trap Index",3, cex=1.5, line=1)
  #   mtext("Trapped per night",3, cex=1.5, line=1)
  #   
  # })
  # 
  
  
  
  output$plot7<-renderPlot({
    par(mar=c(4,3,1,1),tcl=-.2, mgp=c(2,0.5,0))
    locshp<-get.loc.shape(input$sigma.mean, input$sigma.sd)
    hist(rlnorm(10000 ,meanlog=locshp$location, sdlog=locshp$shape), xlab="", ylab="", main="", col="grey")
    mtext("Sigma (m)",1, cex=1, line=2)
    
  })
  
  
  output$plot6<-renderPlot({
    par(mar=c(4,3,1,1),tcl=-.2, mgp=c(2,0.5,0))
    alpbet<-get.alpha.beta(input$g0.mean, input$g0.sd)
    hist(rbeta(10000, alpbet$alpha, alpbet$beta), xlim=c(0,1), main="", col="grey", xlab="", ylab="")
    mtext("g0",1, cex=1, line=2)
    
    
  })
  
  
  output$plot.eff<-renderPlot({
    p.hunt.day<-mydata.zone()$p.hunt.day
    shp.2<-mydata.zone()$shp.2
    eff.id<-cut(p.hunt.day,breaks=seq(from=0, to=1, by=0.125), labels=c(1,2,3,4,5,6,7,8))
    plot(shp.2, main="Prob. kill/day", col=cols.eff[eff.id])
    
    
    legend("topleft",legend=c(paste0("Zone A ", round(p.hunt.day[1],2)), paste0("Zone B ", round(p.hunt.day[2],2)), 
                              paste0("Zone C ", round(p.hunt.day[3],2)), paste0("Zone D ", round(p.hunt.day[4],2))), bty='n')
    
    
    
  })
  
  
  output$plot.hunt<-renderPlot({
    eff<-seq(from=50, to=3000, by=50)
    theta.hat<-1-exp(-((input$hunt.rho*log(eff))^input$hunt.k))
    plot(eff, theta.hat, ylim=c(0,1), ylab="Pr.Kill",las=1, xlab="Distance (m)", type='l')
  })
  
  
  # output$text1<-renderText({
  #   traps<-mydata()$traps
  #   return(paste0("Total traps: ", dim(traps)[1]))
  # })
  # 
  # output$text2<-renderText({
  #   animals<-datab()$animals.xy
  #   return(paste0("Total animals: ", dim(animals)[1]))
  # })
  # 
  # output$text3<-renderText({
  #   trap.catch.mat<-datab()$trap.catch.mat
  #   hunt.catch.mat<-datab()$hunt.catch.mat
  #   a<-rowSums(trap.catch.mat)+rowSums(hunt.catch.mat)
  #   return(paste0("Animals killed: ", mean(a)))
  # })
  
  output$text4<-renderText({
    trap.catch<-datab()$trap.catch
    animals<-datab()$animals.xy
    return(paste0("Prob capture: ", round(sum(trap.catch)/dim(animals)[1],3)))
  })
  
  output$text5<-renderText({
    traps<-mydata.map()$traps
    ntraps<-dim(traps)[1]
    ha.area<-mydata.shp()$ha
    
    paste(ntraps, " Traps\nOne trap per ",round(ha.area/ntraps,1)," hectares", sep="")
  })
  
  
  output$text6<-renderText({
    hr.radius<-input$sigma.mean*2.45
    # ha.area<-input$area.ha
    return(paste0("Home range size: ",round(((pi*hr.radius^2)/10000),2), " ha"))
  })
  
  # output$text7<-renderText({
  #   n.nights<-input$n.nights
  #   pop.size<-datab()$pop.size.mat
  #   return(paste0("Final population size\nMean: ", round(mean(pop.size[,n.nights+1]),2),"\nStdDev: ",round(sd(pop.size[,n.nights+1]),2),
  #                 "\nProbability of Eradication:", round(mean(pop.size[,n.nights+1]==0),2)))
  # })
  
  output$text9<-renderText({
    
    sigma.mean<-input$sigma.mean
    grid.size<-sigma.mean*4
    
    return(paste0("Suggested grid size: ", round(grid.size,0)," m"))
  })
  
  
  output$text8<-renderText({
    
    g0.mean<-input$g0.mean
    max.sd<-sqrt((g0.mean^2 - g0.mean)*-1)
    return(paste0("Max SD: ", round(max.sd,2)))
  })
  
  
  output$text10<-renderText({
    
    ha<-mydata.shp()$ha
    
    return(paste0("Total area (ha): ", round(ha,0)))
  })
  
  # output$text11<-renderText({
  #   ha<-mydata()$ha
  #   hunt.cell.size<-input$hunt.cell.size
  #   return(paste0("Number of cells ", round(ha/hunt.cell.size,0)))
  # })
  
  
  # output$scenario_dropdown<-renderUI({
  #   
  #   params<-mydata.c()$params
  #   
  #   selectInput("result_scenario","Choose Scenario to plot",params$Scenario)
  #   
  # })
  # output$text12<-renderText({
  #   p.hunt.day<-mydata()$p.hunt.day
  #   return(paste0("Daily kill prob\nZone A ", round(p.hunt.day[1],2), "\nZone B ", round(p.hunt.day[2],2), 
  #                 "\nZone C ", round(p.hunt.day[3],2), "\nZone D ", round(p.hunt.day[4],2)))
  # })
  # 
  
  
  # output$txt_trap_cost<-renderText({
  #   a<-mydata()$trap.cost
  #   paste("Trapping cost: $",format(round(a,2), big.mark=",", scientific=FALSE), sep="")
  # })
  # 
  # 
  # output$txt_hunt_cost<-renderText({
  #   day.rate<-input$day.rate.hunt
  #   days.zone<-c(input$hunt.days.a, input$hunt.days.b, input$hunt.days.c, input$hunt.days.d)
  #   a<-day.rate*sum(days.zone)
  #   paste("Hunting cost: $",format(round(a,2), big.mark=",", scientific=FALSE), sep="")
  # })
  # 
  # 
  # 
  # 
  # output$txt_total_cost<-renderText({
  #   a<-mydata()$trap.cost
  #   day.rate<-input$day.rate.hunt
  #   days.zone<-c(input$hunt.days.a, input$hunt.days.b, input$hunt.days.c, input$hunt.days.d)
  #   b<-day.rate*sum(days.zone)
  #   paste("Total cost: $",format(round(a+b,0), big.mark=",", scientific=FALSE), sep="")
  # })
  
  
  
  
  output$results.table<-renderDataTable({
    
    datab()$params
    
  })
  
  output$scenarios.table<-renderTable({
    
    mydata.scen()$params
    
  })
  
  
  #   output$my_table <- DT::renderDataTable(
  #     dat, selection = "none", 
  #     options = list(searching = FALSE, paging=FALSE, ordering=FALSE, dom="t"), 
  #     server = FALSE, escape = FALSE, rownames= FALSE, colnames=c("", ""), 
  #     callback = JS("table.rows().every(function(i, tab, row) {
  #                   var $this = $(this.node());
  #                   $this.attr('id', this.data()[0]);
  #                   $this.addClass('shiny-input-container');
  # });
  #                   Shiny.unbindAll(table.table().node());
  #                   Shiny.bindAll(table.table().node());")
  #   )
  #   
  #   output$test <- renderText({
  #     as.character(input$x11)
  #   })
  
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         The user interface
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui<-fluidPage(theme=shinytheme("flatly"),
              title="Trapping",
              tags$head(
                tags$style(HTML("
                                .shiny-output-error-validation {
                                color: blue;
                                }
                                "))),
              fluidRow(
                column(width=7,  
                       h1("TrapSim: Island Conservation"),
                       "This trapping simulation app is intended to provide guidance as to the approximate amount of trapping/hunting that you may require to achieve a various levels of pest reduction.",
                       # "You can specify an area of a specified size, or upload a shapefile of the area of interest.",
                       # "When you have set up the area, trap layout and pest parameters, click the", strong("Run Trap Sim"), "button at the bottom of the page to run the trapping simulation.",
                       "For background on TrapSim, ", 
                       tags$a(href="https://www.landcareresearch.co.nz/publications/newsletters/kararehe-kino/issue-32/trapsim-an-online-tool-to-help-managers-decide-on-a-trapping-regime", "Click here.", target="_blank"), p(),
                       "For input parameters with red labels, enter values separated by a slash, e.g. 1000/500 ",
                       p(),
                       strong("Hover cursor over each of the input boxes for pop-up help.")
                       
                       # h3("Trapping Simulation Tool")
                ),
                column(4,
                       img(src="manaaki_logo.png", height = 90, align="right", hspace=20,vspace=10),
                       img(src="ari_logo.jpg", height = 90, align="right", hspace=20,vspace=10),
                       img(src="ciss_logo.jpg", height = 90, align="right", hspace=20,vspace=10)
                )
              ),
              
              fixedRow(
                column(width=4,
                       tabsetPanel(
                         tabPanel("1. Area and Pest parameters",
                                  column(width=8,
                                         h4("Area"),
                                         wellPanel(
                                           #   div(style="display:inline-block",
                                           #       tags$div(title="'Specify-Size': specify a total size for a random square; 'Upload Shapefile': upload a shapefile of your study area",
                                           #                
                                           #                radioButtons(inputId="area_type", label="Chose  area type",choices=c("Specify Size (ha)"="Area","Upload Shapefile"="Map"), selected="Map"))),
                                           #   
                                           #   conditionalPanel(
                                           #     condition="input.area_type=='Area'",
                                           #     numericInput(inputId = "area.ha", label="Area (ha)", value=10000, width="120px")
                                           #   ),
                                           #   
                                           #   conditionalPanel(
                                           #     condition="input.area_type=='Map'",
                                           #     
                                           #     tags$div(title="Be sure to select all the components (.shp, .dbf etc)",
                                           #              fileInput(inputId = "shp.file", label="Chose the VCZ Shapefile", accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj"), multiple=TRUE, width="200px")
                                           #     )
                                           #   ),
                                           verbatimTextOutput("text10")
                                         ),
                                         h4("Pest parameters"),
                                         wellPanel(
                                           div(style="display:inline-block;vertical-align:top",
                                               numericInput(inputId = "numb.poss", label="Number (total)", value=100, width="180px")
                                           ),
                                           div(style="display:inline-block;vertical-align:top",
                                               # fileInput(inputId = "ras.1", label="Ascii file of relative abundance", accept=c('.asc'), multiple=FALSE, width="250px")
                                               radioButtons(inputId = "ras.1", label="Relative abundance", choices=c("Random","Habitat Specific"), selected="Habitat Specific",width="250px")
                                           ),
                                           tags$div(title="Sigma x 2.45 is the radius of a circle  where an indivudual spends 95% of its time.",
                                                    h5(strong("Home range (sigma)"))),
                                           div(style="display:inline-block",
                                               numericInput(inputId = "sigma.mean", label="Mean", value=100, width="120px")),
                                           div(style="display:inline-block",
                                               numericInput(inputId="sigma.sd", label='StdDev', value=5, width="120px")),
                                           div(style="display:inline-block;vertical-align:bottom",
                                               checkboxInput(inputId = "show_sigma",label="Show Sigma dist.", value=FALSE)
                                           ),
                                           conditionalPanel(
                                             condition="input.show_sigma==1", 
                                             plotOutput(outputId = "plot7", width = "250px", height="200px")),
                                           verbatimTextOutput("text6"),
                                           
                                           
                                           
                                           h5(strong("Reproductive parameters")),
                                           div(style="display:inline-block",
                                               tags$div(title="This is the maximum rate of increase.",
                                                        numericInput(inputId = "rmax.poss", label="Rmax", value=0.4, width="120px")
                                               )),
                                           div(style="display:inline-block",
                                               tags$div(title="This is the start day of the reproductive period (i.e. start of the birth pulse).",
                                                        numericInput(inputId = "rep.start", label="Start day", value=100, width="120px")
                                               )),
                                           div(style="display:inline-block",
                                               tags$div(title="This is the number of days the reproductive period lasts (i.e. duration of the birth pulse).",
                                                        numericInput(inputId = "rep.nights", label="Length (days)", value=60, width="120px")
                                               ))
                                         )
                                  )
                         ),
                         tabPanel("2. Control Methods",
                                  # column(width=6,
                                  
                                  tabsetPanel(
                                    # column(width=12,
                                    tabPanel("Trapping",
                                             h4("Trapping"),
                                             wellPanel(
                                               tags$div(title="g0 is the probability of trapping an animal on a single night if the trap is at the middle of its home range.",
                                                        h5(strong("Trap probability (g0)"))),
                                               div(style="display:inline-block;vertical-align:bottom",
                                                   numericInput(inputId = "g0.mean", label="Mean", value=0.1, width="120px")),
                                               div(style="display:inline-block;vertical-align:bottom",
                                                   numericInput(inputId="g0.sd", label='StdDev', value=.01, width="120px")),
                                               div(style="display:inline-block;vertical-align:bottom",
                                                   verbatimTextOutput("text8")),
                                               div(style="display:inline-block;vertical-align:bottom",
                                                   checkboxInput(inputId = "show_g0",label="Show g0 dist.", value=FALSE)
                                               ),
                                               conditionalPanel(
                                                 condition="input.show_g0==1", 
                                                 plotOutput(outputId = "plot6", width = "250px", height="200px")),
                                               
                                               radioButtons(inputId="trap_type", label="Chose traps",choices=c("Simulate Traps"="Sim","Upload Locations"="Loc"), 
                                                            selected="Sim"),
                                               conditionalPanel(
                                                 condition="input.trap_type=='Sim'",
                                                 
                                                 div(style="display:inline-block",
                                                     tags$div(id="redtitle",title="The trap spacing in the east-west direction",
                                                              # numericInput(inputId = "traps.x.space", label="Trap space E-W (m)", value="1000",width="135px"))),
                                                              textInput(inputId = "traps.x.space", label="Trap space E-W (m)", value="1000",width="135px"))),
                                                 div(style="display:inline-block",
                                                     tags$div(id="redtitle",title="The trap spacing in the north-south direction",
                                                              textInput(inputId = "traps.y.space", label="Trap space N-S (m)", value="1000",width="135px"))),
                                                 div(style="display:inline-block",
                                                     tags$div(title="The buffer from the edge",
                                                              numericInput(inputId = "traps.buff", label="Edge buffer (m)", value="100", min=0, max=1000, width="110px")))
                                               ),
                                               
                                               conditionalPanel(
                                                 condition="input.trap_type=='Loc'",
                                                 fileInput(inputId = "trap.file", label="Chose the csv with Locations (E-N)", accept=c('.csv'), multiple=FALSE, width="250px")
                                               ),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Nightly probability of by-catch, false triggers etc  ",
                                                            numericInput(inputId = "p.bycatch", label="Daily bycatch rate", value=0, width="130px"))),
                                               
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Maximum catch per trap  ",
                                                            numericInput(inputId = "max.catch", label="Max catch", value=1, width="120px"))),
                                               p(),
                                               div(style="display:inline-block;vertical-align:middle",h5("All Traps")),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(id="redtitle",title="Start night of trapping.",
                                                            textInput(inputId = "trap.start", label="Start night", value="5", width="120px"))),
                                               
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(id="redtitle",title="Number of nights traps are set for.",
                                                            textInput(inputId = "trap.nights", label="Duration (nights)", value="20", width="120px"))),
                                               
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(id="redtitle",title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                            textInput(inputId = "n.check", label="Check interval", value="1", width="120px"))),
                                               
                                               
                                               p(),
                                               div(style="display:inline-block;vertical-align:middle",h5("Zone A")),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Start night of trapping.",
                                                            numericInput(inputId = "trap.start.a", label="Start night", value="5", width="120px"))),

                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Number of nights traps are set for.",
                                                            numericInput(inputId = "trap.nights.a", label="Duration (nights)", value="20", width="120px"))),

                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                            numericInput(inputId = "n.check.a", label="Check interval", value="1", width="120px"))),

                                               p(),
                                               div(style="display:inline-block;vertical-align:top",h5("Zone B")),
                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="Start night of trapping.",
                                                            numericInput(inputId = "trap.start.b", label=NULL, value="5", width="120px"))),

                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="Number of nights traps are set for.",
                                                            numericInput(inputId = "trap.nights.b", label=NULL, value="20", width="120px"))),

                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                            numericInput(inputId = "n.check.b", label=NULL, value="1", width="120px"))),

                                               p(),
                                               div(style="display:inline-block;vertical-align:top",h5("Zone C")),
                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="Start night of trapping.",
                                                            numericInput(inputId = "trap.start.c", label=NULL, value="5", width="120px"))),

                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="Number of nights traps are set for.",
                                                            numericInput(inputId = "trap.nights.c", label=NULL, value="20", width="120px"))),

                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                            numericInput(inputId = "n.check.c", label=NULL, value="1", width="120px"))),

                                               p(),
                                               div(style="display:inline-block;vertical-align:top",h5("Zone D")),
                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="Start night of trapping.",
                                                            numericInput(inputId = "trap.start.d", label=NULL, value="5", width="120px"))),

                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="Number of nights traps are set for.",
                                                            numericInput(inputId = "trap.nights.d", label=NULL, value="20", width="120px"))),

                                               div(style="display:inline-block;vertical-align:top",
                                                   tags$div(title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                            numericInput(inputId = "n.check.d", label=NULL, value="1", width="120px"))),
                                               
                                               tags$style(type="text/css", "#redtitle {color: red}")
                                             )
                                    ),
                                    tabPanel("Hunting",
                                             h4("Hunting"),  
                                             wellPanel(
                                               
                                               div(style="display:inline-block;vertical-align:middle",h5("Zone A")),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.start.a", label="Start day", value="30", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.days.a", label="Days hunted", value="10", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "effort.a", label="Effort per day (m)", value="250", width="150px"))),
                                               p(),
                                               div(style="display:inline-block;vertical-align:middle",h5("Zone B")),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.start.b", label=NULL, value="30", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.days.b", label=NULL, value="10", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "effort.b", label=NULL, value="200", width="150px"))),
                                               p(),
                                               div(style="display:inline-block;vertical-align:middle",h5("Zone C")),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.start.c", label=NULL, value="40", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.days.c", label=NULL, value="10", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "effort.c", label=NULL, value="150", width="150px"))),
                                               p(),
                                               div(style="display:inline-block;vertical-align:middle",h5("Zone D")),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.start.d", label=NULL, value="50", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.days.d", label=NULL, value="30", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "effort.d", label=NULL, value="100", width="150px"))),
                                               
                                               p(),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.k", label="K", value="3.2", width="150px"))),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "hunt.rho", label="rho", value="0.1", width="150px"))),
                                               
                                               plotOutput(outputId = "plot.eff", width = "450px", height="400px"),
                                               plotOutput(outputId = "plot.hunt", width = "450px", height="400px")
                                             )
                                    ),
                                    tabPanel("Costs",
                                             wellPanel(
                                               h4("Trapping"),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Labour cost - day rate ($)",
                                                            numericInput(inputId = "day.rate", label="Trapping day rate ($)", value=400, width="150px")
                                                   )),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Number of traps checked per day",
                                                            numericInput(inputId = "traps.per.day", label="Traps checked per day", value=40, width="150px")
                                                   )),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="Fixed cost ($) per trap",
                                                            numericInput(inputId = "cost.per.trap", label="Fixed cost per trap ($)", value=20, width="150px")
                                                   )),
                                               h4("Hunting"),
                                               div(style="display:inline-block;vertical-align:middle",
                                                   tags$div(title="help text ",
                                                            numericInput(inputId = "day.rate.hunt", label="Hunting day rate ($)", value=500, width="150px")))
                                               
                                               
                                             )
                                    )  #End of column
                                  )
                         ),
                         tabPanel("3. Run and Results",
                                  column(width=4,
                                         wellPanel(
                                           
                                           radioButtons(inputId="sim_type", label="Simulation Type",choices=c("Individual traps"="individ", "Trap density per grid"="grid"), selected="individ", width="200px"),
                                           conditionalPanel(
                                             condition="input.sim_type=='grid'",
                                             verbatimTextOutput("text9"),
                                             numericInput(inputId = "cell.width", label="Grid cell width (m)", value=500, width="150px")
                                           ),#))                         ),
                                           
                                           
                                           div(style="display:inline-block;vertical-align:middle",
                                               numericInput(inputId = "n.nights",label="Simulation length", value=100, width="120px")),
                                           div(style="display:inline-block;vertical-align:middle",
                                               numericInput(inputId = "n.its",label="Iterations", value=5, width="120px")),
                                           div(style="display:inline-block;vertical-align:middle",
                                               actionButton("act.btn.trapsim","Run TrapSim"))
                                           
                                           # ),
                                           # fluidRow(
                                           #   textOutput("text1"),
                                           #   textOutput("text2"),
                                           #   textOutput("text3"),
                                           #   verbatimTextOutput("text7")
                                         )
                                  ),
                                  column(width=8,
                                         fluidRow(
                                           plotOutput(outputId = "plot2", width = "550px", height="350px")  #Cumulative Captures
                                         ),
                                         fluidRow(
                                           plotOutput(outputId = "plot2.hunt", width = "550px", height="350px")
                                         ),
                                         # fluidRow(
                                         #   plotOutput(outputId = "plot5", width = "450px", height="350px")
                                         # ),
                                         
                                         # ),
                                         # column(width=4,
                                         fluidRow(
                                           plotOutput(outputId = "plot4", width = "550px", height="350px") #Population Size
                                         ),
                                         fluidRow(
                                           plotOutput(outputId = "plotzone", width = "550px", height="350px") #Population Size
                                         )
                                         
                                  )
                                  
                         )
                         
                       )),
                column(width=8,
                       tabsetPanel(
                         tabPanel("1. Map of area, traps and animals",
                                  "This is provided to visualise the area. The parameter values specified here are ", strong("not"), " included in the simulations",
                                  p(),
                                  div(style="display:inline-block",
                                      tags$div(title="The trap spacing in the east-west direction",
                                               numericInput(inputId = "traps.x.space.i", label="Trap space E-W (m)", value="1000",width="135px"))),
                                  div(style="display:inline-block",
                                      tags$div(title="The trap spacing in the north-south direction",
                                               numericInput(inputId = "traps.y.space.i", label="Trap space N-S (m)", value="1000",width="135px"))),
                                  div(style="display:inline-block",
                                      tags$div(title="The buffer from the edge",
                                               numericInput(inputId = "traps.buff.i", label="Buffer (m)", value="100", min=0, max=1000, width="110px"))),
                                  div(style="width:300px;display:inline-block",
                                      verbatimTextOutput("text5")),
                                  p(),
                                  div(style="display:inline-block;vertical-align:top",
                                      numericInput(inputId = "numb.poss.i", label="Population size", value=10, width="180px")
                                  ),
                                  
                                  
                                  leafletOutput(outputId = "mymap", height = "1000px", width="1400px")
                         ),
                         tabPanel("2. Scenarios ",
                                  tableOutput('scenarios.table')
                                  
                         ),
                         tabPanel("3. Results ",
                                  # tableOutput('results.table')
                                  dataTableOutput('results.table'),
                                  # uiOutput("scenario_dropdown")
                                  selectInput("result_scenario","Choose Scenario to plot",choices=NULL)
                                  
                         )
                         
                         
                         
                       )
                       
                       
                )
                
                
                
                
              ), #End of Row
              
              # img(src="logolandcare.png", height = 40, align="centre", hspace=20),
              # img(src="MWlogo.jpg", height = 100, align="centre", hspace=20,vspace=10),
              
              h6("v1.5.2: 26 November 2020"),
              h6("Developed by Manaaki Whenua - Landcare Research."),
              h6("email: gormleya@landcareresearch.co.nz")
)


shinyApp(ui=ui, server=server)

#1.2 includes ability to read in a csv file of traps...
#1.3 - multi-capture traps...
#1.4 - costs and rearrange interface...
#1.5 - run scenarios across some parameter combinations
#1.5.2 - outputs by zones as well. 








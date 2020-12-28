#################################################
#               Packages                        #
#################################################

library(MASS)
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(PPtreeViz)
library(gridExtra)
library(reshape2)
library(PPforest)
library(plyr)
library(dplyr)
library(purrr)
library(xtable)
library(tidyr)
library(rpart)
library(microbenchmark)
library(tidyr)
library(randomForest)
library(forcats)
library(patchwork)
library(viridis)

#################################################
#      grid_arrange_shared_legend function      #
#################################################

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.draw(combined)
}

#####################################################
# Simulated data to plot CART and PPtree boundaries #
#####################################################

simu3 <- function(mux1, mux2, muy1, muy2, muz1, muz2,
                  cor1, cor2, cor3) {
  bivn <- mvrnorm(100, mu = c(mux1, mux2),
                  Sigma = matrix(c(1, cor1, cor1, 1), 2))
  bivn2 <- mvrnorm(100, mu = c(muy1, muy2),
                   Sigma = matrix(c(1, cor2, cor2, 1), 2))
  bivn3 <- mvrnorm(100, mu = c(muz1, muz2),
                   Sigma = matrix(c(1, cor3, cor3, 1), 2))
  
  d1 <- data.frame(Sim = "sim1", bivn, stringsAsFactors = TRUE)
  d2 <- data.frame(Sim = "sim2", bivn2, stringsAsFactors = TRUE)
  d3 <- data.frame(Sim = "sim3", bivn3, stringsAsFactors = TRUE)
  rbind(d1, d2, d3)
}

set.seed(666)
dat.pl2 <- simu3(-1,0.6,0,-0.6, 2,-1,0.95, 0.95, 0.95)
dat.plnew <- simu3(0,0.6,0,-0.6, 0,-1,-0.1, -0.1, -0.1)

##################################################
#        Figure 1,  Huber plot                   #
##################################################
Huberplot2 <- function (origdata2D, origclass, PPmethod = "LDA", weight = TRUE, 
                        r = 1, lambda = 0.5, opt.proj = TRUE, UserDefFtn = NULL, stdd = TRUE,
                        ...) 
{
  index <- NULL
  best.proj <- NULL
  best.index <- 0
  origdata2D <- as.matrix(origdata2D)
  for (i in 0:360) {
    theta <- pi/180 * i
    proj.data <- matrix(cos(theta) * origdata2D[, 1] + sin(theta) * 
                          origdata2D[, 2])
    proj <- matrix(c(cos(theta), sin(theta)), ncol = 1)
    if (PPmethod == "LDA") {
      newindex <- LDAindex(origclass, origdata2D, proj = proj, 
                           weight = weight)
    }
    else if (PPmethod == "PDA") {
      newindex <- PDAindex(origclass, origdata2D, proj, 
                           weight = weight, lambda = lambda)
    }
    else if (PPmethod == "Lr") {
      newindex <- Lrindex(origclass, origdata2D, proj, 
                          weight = weight, r = r)
    }
    else if (PPmethod == "GINI") {
      newindex <- GINIindex1D(origclass, origdata2D, proj)
    }
    else if (PPmethod == "ENTROPY") {
      newindex <- ENTROPYindex1D(origclass, origdata2D, 
                                 proj)
    }
    else if (PPmethod == "UserDef") {
      newindex <- UserDefFtn(proj.data, ...)
    }
    index <- c(index, newindex)
  }
  sel.index <- which(index[1:360] > signif(max(index), 6) - 
                       1e-06)
  theta.best.all <- pi/180 * (sel.index - 1)
  theta.best <- theta.best.all[1]
  proj.data.best <- matrix(cos(theta.best) * origdata2D[, 1] + 
                             sin(theta.best) * origdata2D[, 2])
  index.best <- max(index)
  rg <- round(max(index) - min(index), 5)
  PPindex <- (index -  0) * 3 + 2 # Raw index value is assumed to be between 0, 1
  if (stdd) {
    if (rg == 0) {
      PPindex <- rep(4, length(index))
    }
    else {
      PPindex <- (index - min(index))/rg * 2 + 3
    }
  }
  index.median <- median(PPindex)
  data.circle <- NULL
  data.index <- NULL
  for (i in 1:361) {
    theta <- pi/180 * (i - 1)
    data.index <- rbind(data.index, c(PPindex[i] * cos(theta), 
                                      PPindex[i] * sin(theta)))
    data.circle <- rbind(data.circle, c(index.median * cos(theta), 
                                        index.median * sin(theta)))
  }
  maxdiff <- max(c(diff(range(origdata2D[, 1])), diff(range(origdata2D[, 
                                                                       2]))))
  orig.scaled <- apply(origdata2D, 2, function(x) (x - mean(x))/maxdiff * 
                         3.5)
  data.cX <- data.circle[, 1]
  data.cY <- data.circle[, 2]
  data.X <- data.index[, 1]
  data.Y <- data.index[, 2]
  plot.data <- data.frame(data.cX, data.cY, data.X, data.Y)
  x <- orig.scaled[, 1]
  y <- orig.scaled[, 2]
  group <- origclass
  point.data <- data.frame(x, y, group)
  min.X <- min(unlist(plot.data))
  max.X <- max(unlist(plot.data))
  P1 <- ggplot(data = plot.data, aes(x = data.X, y = data.Y)) + 
    geom_path() + 
    geom_path(aes(x = data.cX, y = data.cY), 
              linetype = "dashed") + 
    geom_point(data = point.data, 
               aes(x = x, y = y, color = group, shape = group)) + 
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(breaks = NULL) + 
    xlab("") + ylab("") + 
    coord_fixed() + 
    theme_bw() + 
    scale_shape_manual(values = c(15, 16, 17),labels = c("g1","g2","g3")) +
    scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), 
                        labels = c("g1","g2","g3")) + 
    labs(colour="Class", shape="Class") +
    theme(panel.border = element_blank())
  if (opt.proj) {
    P1 <- P1 + 
      geom_abline(intercept = 0, 
                  slope = sin(theta.best)/cos(theta.best), 
                  linetype="dotted") +
      ggtitle(paste0("Max: ", round(index.best, 2)))
    if (length(theta.best.all) > 1) 
      for (i in 2:length(theta.best.all)) 
        P1 <- P1 + geom_abline(intercept = 0, 
                               slope = sin(theta.best.all[i])/cos(theta.best.all[i]), 
                               linetype = "dotted")
  }
  best.proj.data <- proj.data.best
  group <- origclass
  hist.data <- data.frame(best.proj.data, group)
  P2 <- ggplot(data = hist.data, aes(x = best.proj.data)) + 
    geom_density(aes(group = group, fill = group, colour = group), alpha=0.5) +
    scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), 
                      labels = c("g1","g2","g3")) +
    scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3"), 
                        labels = c("g1","g2","g3")) +
    labs(fill = "Class", colour = "Class", x = "Optimal data projection", y="Count") + theme(legend.position = "none")
  list(P1, P2)
  
}

hu1 <- Huberplot2(dat.plnew[,2:3], dat.plnew[,1], PPmethod = "LDA", stdd = FALSE)


hu2 <- Huberplot2(dat.pl2[,2:3], dat.pl2[,1], PPmethod = "LDA", stdd = FALSE)

grid.arrange(hu1[[1]], hu1[[2]], hu2[[1]], hu2[[2]], ncol = 2)

##################################################
#         Plot CART vs PPtree boundaries         #
##################################################

grilla <- expand.grid(X1 = seq(-4,4.8,,100),
                      X2 = seq(-4.3,3.3,,100))

pptree <- PPtreeViz::PPTreeclass(Sim~., data = dat.pl2, "LDA")
ppred.sim <- PPtreeViz::PPclassify(pptree, test.data = grilla,
                                   Rule = 6)
grilla$ppred<-ppred.sim[[2]]

rpart.crab <- rpart(Sim ~ X1 + X2, data = dat.pl2)
rpart.pred <- predict(rpart.crab, newdata = grilla, type="class")

bnds <- data.frame(
  a = c(pptree$splitCutoff.node[[1]] / pptree$projbest.node[[3]],
        pptree$splitCutoff.node[[2]] / pptree$projbest.node[[4]]),
  b = c(-pptree$projbest.node[[1]] / pptree$projbest.node[[3]],
        -pptree$projbest.node[[2]] / pptree$projbest.node[[4]]))
p <- ggplot(data = grilla ) +
  geom_point(aes(x = X1, y = X2,
                 color = as.factor(ppred),
                 shape=as.factor(ppred)),
             alpha = .20) +
  scale_colour_brewer(name="Class", type="qual",
                      palette = "Dark2") +
  
  scale_shape_discrete(name = 'Class') +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position="none")

pl.pp <- p + geom_point(data = dat.pl2,
                        aes(x = X1, y = X2,
                            group = Sim,
                            shape = Sim, color=Sim), size = 3) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))

p2 <- ggplot(data = grilla) +
  geom_point(aes(x = X1, y = X2,
                 color = as.factor(rpart.pred),
                 shape =  as.factor(rpart.pred)), alpha = .2) +
  scale_colour_brewer(name = "Class",
                      labels = levels(dat.pl2$Sim),
                      type = "qual",palette = "Dark2") +
  theme_bw() + scale_shape_discrete(name = 'Class') +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme(aspect.ratio = 1, legend.position ="none")

pl.rpart <- p2 + geom_point(data = dat.pl2,
                            aes(x = X1 , y = X2,
                                group = Sim, shape = Sim,
                                color = Sim), size = 3)

grid.arrange(pl.rpart, pl.pp, ncol = 2)


##################################################
# To generate file preformance_timesWTplyrnew.Rdata #
# Study PPforest implementation time             #
##################################################

# Functions

library(plyr)

dt_fn <- function(g, n, p) {
  # g: groups; n: obs per group; p: predictor variables,
  x <- matrix(rnorm(p*n*g), ncol = p, nrow = n*g)
  data.frame(type = as.factor( rep(1:g, each = n) ), x)
}

bsq_fn <- function( dd, mtree, cr, prop.v) {
  # dd: data, mtree: number of trees, #cr number of parallel cores, prop.v: prop of variables used
  p <- ncol(dd[,-1])
  PPforest(data = dd, class = "type",
           size.tr = 1, m = mtree , size.p = prop.v,
           PPmethod = 'PDA')
}



micro_fn <- function(gs, ns, ps, prop.vs,mtrees, crs,  dd.list) {
  # identify dataset
  ids <- attributes(dd.list)$split_labels %>%
    mutate(id = 1:length(dd.list)) %>%
    filter(g == gs, n == ns, p == ps)
  
  # run forest
 
    md = microbenchmark( bsq_fn(dd = dts[[ ids$id ]], mtree = mtrees, prop.v = prop.vs), times = 5)
  
  return(md)
}


# datasets


prs.dt <- expand.grid(g = c(3, 6, 9), n = c(10, 100), p = c(10, 100))
dts <- mlply(prs.dt, dt_fn)

# set up scenarios, separating old and new
prs.new <- expand.grid(gs = c(3, 6, 9), ns = c(10, 100),
                       ps = c(10, 100), mtrees = c(50, 500), prop.vs = c(.2, .4, .8) )

pt <- proc.time()
# run NEW scenarios
library(microbenchmark)
library(PPforest)
mm.new <- mdply(prs.new, micro_fn, dd.list = dts)

pt <- proc.time() - pt

save(mm.new, file = 'preformance_timesWTplyrnew.Rdata')

###################################################
#  Time performance PPforest versions plot        #
###################################################
load("preformance_timesWTplyrnew.Rdata")

global_labeller <- labeller(
  g = class,
  ntrees = trees,
  .default = label_both
)

mm.allbig <- mm.new %>%
  mutate(seconds = time/1e9, prop.vs =as.factor(prop.vs))

mm.allbig <- mm.allbig %>%
  rename(n = ns, p = ps, g = gs, B=mtrees)
ggplot(data = mm.allbig) +
  geom_smooth(aes(g, seconds, color = prop.vs,
                  linetype = prop.vs), method = "lm") +
  geom_jitter(aes(g, seconds, color = prop.vs, shape = prop.vs),
              alpha = .4, size = 3, height = 0) +
  scale_x_continuous(breaks=seq(0,10,2)) +
  labs(x = "Num. groups", y = "Time (sec)") +
  scale_colour_viridis_d(name = "Proportion of variables", 
                         begin = 0.1, end = 0.6, 
                         option = "magma") +
  scale_shape(name = "Proportion of variables") +
  scale_linetype(name = "Proportion of variables") +
  facet_grid(B~n + p, labeller = label_both,
             scales = "free_y") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 6), aspect.ratio = 1 )

#################################################
#            Summary of benchmark data          #
#################################################
d <- data(package = "PPforest")
# names of the datasets: d$results[, "Item"]
nm <- d$results[, "Item"]
## call the promised data
data(list = nm, package = "PPforest")
## get the dimensions of each data set
tbl2 <- map_df(mget(nm), function(dd) {
  ng <- as.numeric(table(dd[,1]))/nrow(dd) # group size weights
  
  cnd <- sapply(dd, is.factor)
  cnd[1] <- TRUE
  R <- cor(dd[ , !cnd])
  tibble(Cases = nrow(dd), Predictors = as.integer(ncol(dd)-1), Groups = length(ng),
         Imbalance = max(ng) - min(ng), Correlation =
           mean( abs(R[lower.tri(R)])) )}) %>%
  mutate(Data = nm) %>% arrange(desc(Correlation)) %>%
  dplyr::select(Data, Cases:Correlation) 

xtable(tbl2, caption = 'Overview of benchmark data: number of cases, predictors, groups, imbalance and correlation. Imbalance indicates relative class sizes (0=balanced classes), and higher correlation indicates the potential for separations occurring in combinations of variables.', label = 'bench.tab') %>%
  print(caption.placement = 'top', include.rownames = FALSE)

##############################################
#           Benchmark data performance      #
#############################################

# to replicate the .Rdata with the benchmarck performance

train_fn <- function(data, class, size.p) {
  id <- NULL
  n <- nrow(data)
  class.id <- data %>%
    dplyr::select_(class) %>%
    dplyr::mutate(id = 1:n)
  
  class.id %>%
    dplyr::group_by_(class) %>%
    dplyr::sample_frac(size.p) %>%
    dplyr::arrange(id) %>%
    dplyr::ungroup() %>%
    dplyr::select(id)
}

#PPtree
pptr.err <- function(tr, te, met){
  ppt <- PPTreeclass(Type~., data = tr, PPmethod = met)
  m.tr <- PPclassify(ppt, test.data = tr[,-1], true.class = tr[,1],  Rule = 1)
  m.te <- PPclassify(ppt, test.data = te[,-1], true.class = te[,1] , Rule = 1)
  data_frame(err.tr = m.tr[[1]]/length(m.tr[[2]]), err.te = m.te[[1]]/length(m.te[[2]]))
  
}

#CART
cart.err <- function(tr, te){
  cart.m <- rpart(as.factor(Type)~., data = tr)
  m.tr <- predict(cart.m, newdata = tr, type = "class")
  m.te <- predict(cart.m, newdata = te[,-1], type = "class")
  tab.tr <- table(tr[ , 1], m.tr)
  tab.te <- table(te[ , 1], m.te)
  data.frame(err.tr = (dim(tr)[1] - sum(diag( tab.tr )))/dim(tr)[1] ,
             err.te = (dim(te)[1] - sum(diag( tab.te )))/dim(te)[1] )
}



#RF
rf.err <- function(tr, te) {
  rf <- randomForest(as.factor(Type) ~ ., data = tr)
  m.te <- predict(rf, newdata = te[,-1], type = "class")
  tab.te <- table(te[,1], m.te)
  data.frame(err.tr = ( dim(tr)[1] - sum(diag(rf$confusion[ , -dim(rf$confusion)[2]]) ) )/dim(tr)[1]
             , err.te =  ( dim(te)[1] - sum(diag(tab.te)))/dim(te)[1])
}


#PP.RF
pprf.err <- function(tr, te, lam = 0.1, met, sp) {
  
  ppf.m <- PPforest(tr, "Type", std = FALSE, size.tr = 1, m = 500, PPmethod = met, size.p = sp,
                    lambda = lam)
  tab.te <- table(te[,1],trees_pred(ppf.m, xnew = te[,-1])[[2]] )
  data.frame( err.tr =  ppf.m$training.error,
              err.te = ( dim(te)[1] - sum(diag(tab.te)))/dim(te)[1] )
}


# d <- leukemia
# r <-5
set.seed(123)
models.errpda <- function(d, r = 200){
  
  d %>% crossing(rep = 1:r) %>%
    group_by(rep) %>% nest() %>%
    mutate(errors = map(data, .f = function(ddt) {
      tr.index <- train_fn(ddt, class = "Type", size.p = 2/3)
      te1 <- ddt[-tr.index$id, ] %>% as.data.frame()
      tr1 <- ddt[tr.index$id, ] %>% as.data.frame()
      list(pptr.pda = pptr.err(tr = tr1, te = te1, met = "PDA"),
           pptr.lda = pptr.err(tr = tr1, te = te1, met = "LDA"),
           cart = cart.err(tr = tr1, te = te1),
           rf = rf.err(tr = tr1, te = te1),
        ppf.pdagen = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = sqrt(ncol(tr1)-1)/(ncol(tr1)-1)),
        ppf.pda.6 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp =.6),
        ppf.pda.9 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp =.9))%>%
         bind_rows( .id = "method" )
    })) %>% dplyr::select(-data) %>% unnest(errors) %>%
    group_by(method) %>% 
    dplyr::summarise( mn.tr = mean(err.tr), mn.te = mean(err.te),
                                    md.tr = median(err.tr), md.te = median(err.te) )
  
}

models.errlda <- function(d, r = 10){
  d %>% crossing(rep = 1:r) %>%
    group_by(rep) %>% nest() %>%
    mutate(errors = map(data, .f = function(ddt) {
      tr.index <- train_fn(ddt, class = "Type", size.p = 2/3)
      te1 <- ddt[-tr.index$id, ] %>% as.data.frame()
      tr1 <- ddt[tr.index$id, ] %>% as.data.frame()
      list(pptr.pda = pptr.err(tr = tr1, te = te1, met = "PDA"),
           pptr.lda = pptr.err(tr = tr1, te = te1, met = "LDA"),
           cart = cart.err(tr = tr1, te = te1),
           rf = rf.err(tr = tr1, te = te1),
           #ppf.ldagen = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = sqrt(ncol(tr1)-1)/(ncol(tr1)-1)),
           ppf.lda.6 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .6),
           ppf.lda.9 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp =.9)) %>%
        bind_rows( .id = "method" )
    })) %>% dplyr::select(-data) %>% unnest(errors) %>%
    group_by(method) %>% 
    dplyr::summarise( mn.tr = mean(err.tr), mn.te = mean(err.te),
                                    md.tr = median(err.tr), md.te = median(err.te) )
  
}


lapply(list.files(pattern='result.RData'), load)
lapply(list.files(pattern='_result.Rdata'), load)

l1 <- list(lymphoma_result, NCI60_result, parkinson_result, fishcatch_result, leukemia_result, olive_result, wine_result,
           image_result, glass_result)


table.raw <- crab_result %>% bind_rows(l1)

#generate the table to get the benchmark results

save(table.raw, file = "table_raw.Rdata")

##############################################
#   Table with benchmark results             #
##############################################

load("table_raw.Rdata")

tt1 <- table.raw %>%
  separate(method, into = c('met','b')) %>%
  group_by(data, met) %>% mutate(mm = min(mn.te)) %>%
  filter(mn.te == mm) %>%
  select(data, met, mn.te) %>% spread(met, mn.te) %>% 
  rename(Data =  data, TEcart = cart, TEppf = ppf, TEpptr = pptr, TErf = rf) %>% 
  select(Data,TEcart, TEpptr, TErf, TEppf) 

#Run this code to  replicates the table below
tt2 <-  table.raw %>%
  separate(method, into = c('met','b')) %>%
  group_by(data, met) %>% mutate(mm = min(mn.te)) %>%
  filter(mn.te == mm) %>%
  select(data, met, mn.tr) %>% spread(met, mn.tr) %>% 
  rename(Data =  data)  %>%
  inner_join(select(tbl2, Correlation, Data) ) %>%
  arrange(desc(Correlation)) %>% 
  select(Data,cart, pptr, rf, ppf) 

tbl3 <- tt2 %>%  inner_join(tt1) %>% 
  rename(CART = cart,  PPtree = pptr, RF = rf, PPforest = ppf,
         CART1 = TEcart, PPtree1 = TEpptr, RF1 = TErf,PPforest1 = TEppf) 


# para los encabezados de columnas
addtorow <- list()
addtorow$pos <- list(-1) # para que quede arriba de los column names
addtorow$command <- c('\\hline &  \\multicolumn{4}{c|}{TRAINING} & \\multicolumn{4}{c}{TEST} \\\\')

tbl3 %>% xtable(align = "ll|cccc|cccc", caption = 'Comparison of  CART, PPtree, RF and PPF results with various data sets. The mean of training and test error rates from 200 re-samples is shown. (Order of rows is same as in Table \\ref{bench.tab}.) PPF performs favorably compared to the other methods. \\label{res}',  digits = 3) %>%
  print(add.to.row = addtorow, include.rownames = F, caption.placement = 'top', 
        sanitize.colnames.function = function(x) gsub("1","",x)) 



#################################################
#             Plot benchmark results            #
#################################################

load("table.raw.Rdata")
tbl1 <- table.raw %>%rename(Data=data) %>% 
  inner_join(select(tbl2, Correlation, Data) ) %>%
  arrange(desc(Correlation)) %>% 
  select(-Correlation) %>% 
  separate(method, into = c('met', 'b') ) %>% 
  group_by(Data, met) %>% mutate(mm = min(mn.te) ) %>%
  filter(mn.te == mm) %>% select(Data, met, mn.tr, mn.te) %>% ungroup() %>%
  gather(type, error, mn.tr, mn.te) %>% mutate(met = reorder(met, error) ) %>%
  unite(gg, Data, type, remove = FALSE) %>% 
  mutate(type = factor(type, labels = c("Test","Training")),
         met = fct_recode(met, "PPforest" = "ppf", "CART" = "cart",
                          "PPtree" = "pptr", "RF" = "rf"), Data = factor(Data, levels = tbl2$Data))  

setFactorOrder <- function(x, order = sort(levels(x))) {
  # Returns a factor ordered by `order`.
  # If order is missing, defaults to `levels(x)` if available, else to `sort(unique(x))`
  # Useful for ggplot and elsewhere were ordering is based on the order of the levels
  
  if (!is.factor(x)) {
    warning("`x` is not a factor. Will coerce.")
    levs <- sort(unique(x))
    if (missing(order))
      order <- levs
  } else {
    levs <- levels(x)
  }
  
  # any values in order, not in levels(x)
  NotInx <- setdiff(order, levs)
  
  if (length(NotInx)) {
    warning ("Some values not in x:\n", paste(NotInx, collapse=", "))
  }
  
  # levels(x) not explicitly named in order
  Remaining <-  setdiff(levs, order)
  
  order <- c(setdiff(order, NotInx), Remaining)
  
  factor(x, level = order)
}

tbl1$type<- setFactorOrder(tbl1$type, c("Training", "Test"))
ggplot(tbl1) + geom_line(aes(x = fct_relevel(met, "CART", "PPtree", "RF", "PPforest"), y = error, group = gg, color = Data) ) +
  scale_x_discrete( expand = c(0.01, 0.01) ) +
  scale_colour_brewer(palette = "Spectral") +
  
  labs(y ="Average error rate", x = "Method", colour = "Data") +
  theme(axis.text.x = element_text(angle = 90)) +
  
  facet_wrap(~type, ncol = 2, labeller = label_parsed) +
  theme(legend.position = "bottom")  



###################################################
#                 Global importance plot          #
###################################################
set.seed(100)

lymphppf <- PPforest(data = lymphoma, class ="Type", std = TRUE, size.tr = 1, m = 500, size.p = 0.6, PPmethod = 'PDA', lambda=0.1)

lymphrf <- randomForest(Type ~ ., data = lymphoma, importance = TRUE, proximity = TRUE, mtry = 6)

globalimpo <- ppf_global_imp(data = lymphoma, class = 'Type', lymphppf)
averimpo <- ppf_avg_imp(lymphppf, "Type")
imporf2 <- data.frame(Type = as.factor(rownames( lymphrf$importance)), lymphrf$importance) %>%
  select(Type, MeanDecreaseAccuracy) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>% top_n(10)

ppfimpo1 <- ggplot(averimpo[1:10,], aes(x = mean, y = variable)) + geom_point() + theme(aspect.ratio = 1) + labs(x = "Global aver. importance (VI2)", title ="PPforest", y="")

ppfimpo2 <- ggplot(globalimpo[1:10,], aes(x = mean, y = variable)) + geom_point() + theme(aspect.ratio = 1) + labs(x = "Global importance (VI3)", title ="PPforest", y ="")

rfimpo <- ggplot(imporf2, aes(y = fct_reorder(Type, MeanDecreaseAccuracy ), x = MeanDecreaseAccuracy)) + geom_point()  + theme(aspect.ratio = 1) + labs(y = "", x = "Permuted importance (VI1)", title="Random Forest")

grid.arrange(rfimpo,ppfimpo1, ppfimpo2,  ncol=3)

p1 <- ggplot(lymphoma, aes(x = Gene35, y = Gene50, colour = Type)) +
  geom_point() + theme(aspect.ratio  = 1) + ggtitle("PPforest") +
  scale_colour_brewer(type = "qual", palette = "Dark2")
p2 <- ggplot(lymphoma, aes(x = Gene35, y = Gene44, colour = Type)) +
  geom_point() + theme(aspect.ratio = 1) + ggtitle("PPforest") +
  scale_colour_brewer(type = "qual", palette = "Dark2")
p3 <- ggplot(lymphoma, aes(x = Gene34, y = Gene49, colour = Type)) +
  geom_point() + theme(aspect.ratio = 1) + ggtitle("Random Forest") +
  scale_colour_brewer(type = "qual", palette = "Dark2")

grid_arrange_shared_legend(p3,p1, p2)

###############################################
# Ternary and side-by-side plot               #
###############################################

f.helmert <- function(d)
{
  helmert <- rep(1 / sqrt(d), d)
  for (i in 1:(d - 1))
  {
    x <- rep(1 / sqrt(i * (i + 1)), i)
    x <- c(x, -i / sqrt(i * (i + 1)))
    x <- c(x, rep(0, d - i - 1))
    helmert <- rbind(helmert, x)
  }
  
  return(helmert)
}

#ppf PPforest object
#V1,V2,V3 select the 3 proj directions
ppf <- lymphppf
rf <- lymphrf
n.class <- ppf$train %>% select_(ppf$class.var) %>% unique() %>% nrow()
projct <- t(f.helmert(length(unique(ppf$train[, ppf$class.var])))[-1,])

datppf <-
  data.frame(
    Class = ppf$train[, ppf$class.var],
    ids = 1:nrow(ppf$train),
    proj.vote = as.matrix(ppf$votes) %*% projct
  )


datrf <-
  data.frame(
    Class = ppf$train[, ppf$class.var],
    ids = 1:nrow(ppf$train),
    proj.vote = as.matrix(rf$votes) %*% projct
  )

#need to run f_hermite
f_composition <- function(data) {
  d <- dim(data)[2]
  hm <- f.helmert(d)
  x <- data - matrix(1 / d, dim(data)[1], d)
  return((x %*% t(hm))[,-1])
}

simplex <- function(p = 3) {
  vert <- f_composition(diag(p + 1))
  colnames(vert) <- paste0("d", 1:ncol(vert))
  
  wires <-
    do.call(expand.grid, list(c(1:nrow(vert)), c(1:nrow(vert))))
  
  structure(list(points = vert,
                 edges = wires[!(wires[, 1] == wires[, 2]),]))
}

##ternary plot
s <- simplex(2)
pts <- data.frame(s$points)

edg <- data.frame(x1=pts[,"d1"][s$edges[,1]], x2=pts[,"d1"][s$edg[,2]],
                  y1=pts[,"d2"][s$edg[,1]], y2=pts[,"d2"][s$edg[,2]])

ternaryppf <- datppf %>% ggplot(aes(proj.vote.x, proj.vote.x.1, color = Class)) +
  geom_segment(data = edg, aes(x = x1, xend = x2,
                               y = y1, yend = y2), color = "black" ) +
  geom_point(size = I(3), alpha = 2/3) +
  labs(y = "",  x = "") +
  theme(legend.position = "bottom", aspect.ratio = 1) +
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  labs(x = "T1", y = "T2", title = "PPforest") +
  theme(aspect.ratio=1)  + scale_y_reverse()

ternaryrf <- datrf %>% ggplot(aes(proj.vote.x, proj.vote.x.1, color = Class)) +
  geom_segment(data = edg, aes(x = x1, xend = x2,
                               y = y1, yend = y2), color = "black" ) +
  geom_point(size = I(3), alpha = 2/3) +
  labs(y = "",  x = "") +
  theme(legend.position = "bottom", aspect.ratio = 1) +
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  labs(x = "T1", y = "T2", title ="Random Forest") +
  theme(aspect.ratio=1)  + scale_y_reverse()



voteinfppf <- data.frame(ids = 1:length(lymphppf$train[, 1]), Type = lymphppf$train[, 1],
                         lymphppf$votes, pred = lymphppf$prediction.oob ) %>%
  gather(Class, Probability, -pred, -ids, -Type)

voteinfrf <- data.frame(ids = 1:length(lymphrf$predicted), Type = lymphrf$y,
                        lymphrf$votes, pred = lymphrf$predicted ) %>%
  gather(Class, Probability, -pred, -ids, -Type)

side <-  function(ppf, voteinf, ang = 0, lege = "bottom", siz = 6,
                  ttl = "Side by side dotplot") {
  
  ggplot(data = voteinf, aes(x = fct_reorder(Class, Probability,,.desc = TRUE), y = Probability, color = Type)) +
    geom_jitter(height = 0, width=0.2, size = I(siz), alpha=2/3) +
    ggtitle(ttl) +
    scale_colour_brewer(type = "qual", palette = "Dark2" ) +
    theme(legend.position = lege, legend.text = element_text(angle = ang), axis.text.x = element_text(angle = 45), aspect.ratio = 1) +
    labs(colour = "Class", y = "Proportion", x = "")
}

sidepp <- side(oliveppf, voteinfppf, ttl = "PPforest", siz = 2, lege = "none")
siderf <- side(oliverf, voteinfrf, ttl = "Random Forest", siz = 2)

grid_arrange_shared_legend( ternaryrf,ternaryppf,  siderf,sidepp, ncol = 2, nrow = 2)
################################################
# MDS with proximity matrix info               #
################################################

k = 2
lege = "none"
siz = 3
d <- diag(nrow(lymphppf$train))
dppf <- as.dist(d + 1 - lymphppf$proximity)
pprf.mds <- stats::cmdscale(dppf, eig = TRUE,  k = k)
colnames(pprf.mds$points) <- paste("MDS", 1:k, sep = "")
dfppf <- data.frame(Class = lymphppf$train[, 1], pprf.mds$points)
p12ppf <-  ggplot2::ggplot(data = dfppf) +
  geom_point(ggplot2::aes(x = MDS1, y = MDS2, color = Class),
             size = I(siz), alpha = .5) +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Class") +
  theme(legend.position = lege, aspect.ratio = 1) + labs(title = "PPforest")

drf <- as.dist(d + 1 - lymphrf$proximity)
rf.mds <- stats::cmdscale(drf, eig = TRUE,  k = k)
colnames(rf.mds$points) <- paste("MDS", 1:k, sep = "")
dfrf <- data.frame(Class = lymphppf$train[, 1], rf.mds$points)

p12rf <-  ggplot2::ggplot(data = dfrf) +
  geom_point(ggplot2::aes(x = MDS1, y = MDS2, color = Class),
             size = I(siz), alpha = .5) +
  scale_colour_brewer(type = "qual", palette = "Dark2", name = "Class") +
  theme(legend.position = lege, aspect.ratio = 1) + labs(title ="Random Forest")

grid_arrange_shared_legend(p12rf,p12ppf, nrow=1)


##############################################
# Optimize number of variables, obtaining    #
# the  mean of training and test error rates #
# from 200 re-samples.                       #
#########################################

set.seed(123)

models.errpdavar <- function(d, r = 20){
  d %>% crossing(rep = 1:r) %>%
    group_by(rep) %>% nest() %>%
    mutate(errors = map(data, .f = function(ddt) {
      tr.index <- train_fn(ddt, class = "Type", size.p = 2/3)
      te1 <- ddt[-tr.index$id, ] %>% as.data.frame()
      tr1 <- ddt[tr.index$id, ] %>% as.data.frame()
      list( ppf.pda.1 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .1),
            ppf.pda.2 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .2),
            ppf.pda.3 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .3),
            ppf.pda.4 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .4),
            ppf.pda.5 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .5),
            ppf.pda.6 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .6),
            ppf.pda.7 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp = .7),
            ppf.pda.8 = pprf.err(tr =  tr1, te = te1, lam = 0.1, met = "PDA", sp = .8),
            ppf.pda.9 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "PDA", sp =.9))%>%
        bind_rows( .id = "method" )
    })) %>% dplyr::select(-data) %>% unnest(errors) %>%
    group_by(method) %>% summarise( mn.tr = mean(err.tr), sdmn.tr = sd(err.tr),
                                    mn.te = mean(err.te), sdmn.te = sd(err.te),
                                    md.tr = median(err.tr), sdmd.tr = sd(err.tr),
                                    md.te = median(err.te), sdmd.te =sd(err.te))
  
}


models.errldavar <- function(d, r = 200){
  d %>% crossing(rep = 1:r) %>%
    group_by(rep) %>% nest() %>%
    mutate(errors = map(data, .f = function(ddt) {
      tr.index <- train_fn(ddt, class = "Type", size.p = 2/3)
      te1 <- ddt[-tr.index$id, ] %>% as.data.frame()
      tr1 <- ddt[tr.index$id, ] %>% as.data.frame()
      list(
        ppf.lda.1 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .1),
        ppf.lda.2 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .2),
        ppf.lda.3 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .3),
        ppf.lda.4 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .4),
        ppf.lda.5 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .5),
        ppf.lda.6 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .6),
        ppf.lda.7 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .7),
        ppf.lda.8 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp = .8),
        ppf.lda.9 = pprf.err(tr = tr1, te = te1, lam = 0.1, met = "LDA", sp =.9)) %>%
        bind_rows( .id = "method" )
    })) %>% dplyr::select(-data) %>% unnest(errors) %>%
    group_by(method) %>% summarise( mn.tr = mean(err.tr), sdmn.tr = sd(err.tr),
                                    mn.te = mean(err.te), sdmn.te = sd(err.te),
                                    md.tr = median(err.tr), sdmd.tr = sd(err.tr),
                                    md.te = median(err.te), sdmd.te =sd(err.te) )
  
}



#LDA ok
crab_resvar <- models.errlda(crab)
crab_resvar <- data.frame(Data = "crab",crab_resvar )
save(crab_resvar , file = "crab_resvar.Rdata")

# PDA ok
leukemia_resultvar <- models.errpdavar(leukemia)
leukemia_resultvar <- data.frame(Data = "leukemia", leukemia_resultvar)
save(leukemia_resultvar , file ="leukemia_resultvar.Rdata")

# PDA ok
lymphoma_resultvar <-  models.errpdavar(lymphoma)
lymphoma_resultvar <- data.frame(Data = "lymphoma", lymphoma_resultvar)
save(lymphoma_resultvar , file = " lymphoma_resultvar.Rdata")

# LDA ok
wine_resultvar <-  models.errldavar(wine)
wine_resultvar <- data.frame(Data = "wine", wine_resultvar)
save(wine_resultvar , file ="wine_resultvar.Rdata")

# PDA ok
glass_resultvar <- models.errpdavar(glass)
glass_resultvar <- data.frame(Data = "glass", glass_resultvar)
save(glass_resultvar , file ="glass_resultvar.Rdata")

# LDA  ok
fishcatch_resultvar <- models.errldavar(fishcatch)
fishcatch_resultvar <- data.frame(Data = "fishcatch", fishcatch_resultvar)
save(fishcatch_resultvar , file = "fishcatch_resultvar.Rata")

# LDA ok
parkinson_resultvar <- models.errldavar(parkinson)
parkinson_resultvar <- data.frame(Data = "parkinson", parkinson_resultvar)
save(parkinson_resultvar , file = "parkinson_resultvar.Rdata")

# PDA ok
image_resultvar <- models.errpdavar(image)
image_resultvar <- data.frame(Data = "image",image_resultvar)
save(image_resultvar, file = "image_resultvar.Rdata")

# PDA ok
olive_resultvar <- models.errpdavar(olive)
olive_resultvar <- data.frame(Data = "olive",  olive_resultvar)
save(olive_resultvar, file="olive_resultvar.Rdata")

#PDA ok
NCI60_resultvar <- models.errpdavar(NCI60)
NCI60_resultvar <- data.frame(Data = "NCI60",NCI60_resultvar)
save(NCI60_resultvar, file = "NCI60_resultvar.Rdata")


examvar <- data.frame(sam.p = rep(seq(0.1, 0.9, 0.1), 10),
                      rbind( crab_resvar, leukemia_resultvar, lymphoma_resultvar, wine_resultvar, glass_resultvar, fishcatch_resultvar, parkinson_resultvar, image_resultvar, olive_resultvar,NCI60_resultvar))
save(examvar, file = "examvar.Rdata")


###############################################
#             plot errorvar                   #
###############################################

load("examvar.Rdata")
examvar %>% ggplot(aes(x =  sam.p, y = mn.te, color = Data)) + 
  geom_line() + geom_point() + theme(legend.position = "bottom") + 
  labs(y= "Mean error", x = "Proportion of variables") +
  scale_color_brewer("", palette = "Paired") -> p1


############################################
#   To get errorrates.csv file            #
############################################

set.seed(100)
#to check the number of trees in the forest for each data set
#run for each data set 

tr <- seq(10, 500, 10)
errorrates <- function(dat, met){
 errppf <- NULL
 errrf <-  NULL
for(i in tr){
 ppf <-   PPforest(data = dat, class = "Type", std = TRUE, size.tr = 1, m = i, size.p = 
                     round(sqrt(ncol(dat) -1)/(ncol(dat) - 1), 2),
                   PPmethod = met, lambda = 0.1)

 rf <- randomForest( Type~., data = dat, importance = TRUE,
                          ntree = i)

 errppf <- c(errppf, ppf$oob.error.forest)
 errrf <- c(errrf, rf$err.rate[i,1])
}

data.frame(Trees = tr, PPFerror = errppf, RFerror = errrf )

}

errrate <-  read.csv("errorrates_lymphoma.csv")
errrate <- errrate %>% rename(PPF=PPFerror, RF=RFerror)
errrate %>% gather(Model,ooberror, -Trees, -X) %>% mutate(Model = factor(Model)) %>%
  ggplot(aes(x = Trees, y = ooberror, color = Model)) +
  geom_point(alpha = 0.6) + geom_smooth(se = FALSE) +
  scale_color_brewer("", palette = "Dark2") +
  ylim(c(0,0.15))+
  theme(legend.position = "bottom") +
  labs(y = "OOB error rate") -> p2
grid.arrange(p1, p2, ncol = 2)


###################################################
#  RNA-seq gene expression aplication              #
###################################################
tab_RNA <- readRDS('resRNA_new.rds')
colnames(tab_RNA) <- c("Method", "Training", "Test")
tbrna<- tab_RNA %>% 
  mutate(Method = fct_recode(Method, "CART" = "cart", 
                             "PPF.4"="ppf.pda.4", "PPF.03"="ppf.pda.03", 
                             "PPtree" = "pptr.pda", "RF"="rf"), 
         Training = Training, Test =  Test) %>% 
  xtable(align = "ll||c||c",caption = 'Comparison of CART, PPF, PPtree, 
         and RF results with maize RNA-seq gene exression data set. The 
         mean of training and test error rates from leave-one-out cross validation
         is shown. PPF.03 performs favorably compared to the other methods. \\label{res_RNA}',  digits =3)

tbrna %>% print( include.rownames=F, caption.placement = 'top' )

rnaseq <- readRDS('rnaseq_data.rds')

dat_rna <- rnaseq %>% 
  dplyr::select(GeneID, genotype, total,replicate) %>%
  spread(GeneID,total) %>% 
  dplyr::select(-replicate) %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  dplyr::rename( Type=genotype  ) %>% data.frame()

dd <- dat_rna
dd_sc <- data.frame(dd$Type, scale(dd[,-1]) )
colnames(dd_sc)[1] <- 'Type'
set.seed(123)

# LEAVE-ONE-OUT
# PPtree
library(PPtreeViz)
pptr.err <- function(idx, dd, met){
  tr <- dd[-idx, ]
  ppt <- PPTreeclass(Type~., data = tr, PPmethod = met)
  pp <- PPclassify(ppt, test.data = dd[,-1], true.class = dd[,1],  Rule = 1)[[2]]
  data.frame(err.tr = mean(pp[-idx]!=tr[,1]),
             err.te = 1*(pp[idx]!=dd[idx,1]) )
  
}

# CART
library(rpart)
cart.err <- function(idx, dd){
  tr <- dd[-idx, ]
  cart.m <- rpart(as.factor(Type)~., data = tr)
  pp <- predict(cart.m, newdata = dd, type = "class")
  data.frame(err.tr = mean(pp[-idx]!=tr[,1]),
             err.te = 1*(pp[idx]!=dd[idx,1]) )
  
}

# RF
rf.err <- function(idx, dd) {
  tr <- dd[-idx,]
  rf <- randomForest(as.factor(Type) ~ ., data = tr)
  pp <- predict(rf, newdata = dd[,-1], type = "class")
  data.frame(err.tr = mean(pp[-idx]!=tr[,1]),
             err.te = 1*(pp[idx]!=dd[idx,1]) )
}

# PP.RF
pprf.err <- function(idx, dd, lam = 0.1, met='PDA', sp=.4) {
  dd_sc <- data.frame(dd$Type, scale(dd[,-1]) )
  colnames(dd_sc)[1] <- 'Type'
  tr <- dd_sc[ -idx, ]
  ppf.m <- PPforest(tr, "Type", std = FALSE, 
                    size.tr = 1, m = 500, PPmethod = met,
                    size.p = sp, lambda = lam)
  
  pp <- trees_pred(ppf.m, xnew = dd_sc[,-1]  )[[2]] %>% as.factor()
  levels(pp) <- levels(dd[,1])
  #table(ppf.m$predicting.training, pp[-idx] )
  
  data.frame(err.tr = mean(pp[-idx]!=tr[,1]),
             err.te = 1*(pp[idx]!=dd[idx,1]) )
}

dd <- dat_rna

ll <- list(
  pptr.pda = t(sapply(1:16, pptr.err, dd=dd, met = "PDA")),
  cart = t(sapply(1:16, cart.err, dd=dd)),
  rf = t(sapply(1:16, rf.err, dd=dd)),
  ppf.pda.4 = t(sapply(1:16, pprf.err, dd=dd, lam = 0.1, met = "PDA", sp =.4)),
  ppf.pda.03 = t(sapply(1:16, pprf.err, dd=dd, lam = 0.1, met = "PDA", sp =.03))
)

res.rna <- lapply(ll, data.frame) %>% bind_rows( .id = 'method') %>% 
  group_by(method) %>% 
  summarise(err.tr = mean(as.numeric(err.tr)), 
            err.te = mean(as.numeric(err.te))
  )

saveRDS(res.rna, 'resRNA_new.rds')  

set.seed(123)
ppf.m <- PPforest(dd_sc, "Type", std = FALSE, 
                  size.tr = 1, m = 500, PPmethod = "PDA",
                  size.p = 0.03, lambda = 0.1)

save(ppf.m,file="ppf_m.Rdata")

globalimpo_rna <- ppf_global_imp(data = dd, class = 'Type', ppf.m)
save(globalimpo_rna,file="globalimpo_rna.Rdata")  
averimpo_rna <- ppf_avg_imp(ppf.m, "Type")
save(averimpo_rna,file="averimpo_rna.Rdata")  

# Permute the response, and check results
dd_sc_p <- dd_sc %>% mutate(Type = sample(Type))
ppf.m_p <- PPforest(dd_sc_p, "Type", std = FALSE, 
                    size.tr = 1, m = 500, PPmethod = "PDA",
                    size.p = 0.03, lambda = 0.1)

globalimpo_rna <- ppf_global_imp(data = dd, class = 'Type', ppf.m)
save(globalimpo_rna,file="globalimpo_rna.Rdata")  
averimpo_rna <- ppf_avg_imp(ppf.m, "Type")
save(averimpo_rna,file="averimpo_rna.Rdata")  

##############################################

#############################################

load("globalimpo_rna.Rdata")
load("averimpo_rna.Rdata")
imp_genes10 <- globalimpo_rna %>% 
  slice(1:10) %>% pull(variable) %>% 
  as.character()


#scaled data
pl_rna1 <- dd_sc  %>%
  dplyr::select(imp_genes10,Type ) %>%
  mutate(id = 1:16 )%>%
  gather(key,value,-id,-Type) %>%
  left_join(globalimpo_rna, by=c('key'='variable' ) ) %>% 
  ggplot( aes(x = fct_reorder(key,-mean), y = value,
              group = id, key = id, colour = Type, var = key))  +
  geom_line() + labs(x = "Most important genes", y="Standardize gene expression") +
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  theme( axis.text.x  = element_blank() ) +
  scale_x_discrete( expand = c(0.01, 0.01)) 


sca_rna1 <- dd_sc  %>%
  dplyr::select(imp_genes10,Type ) %>%
  ggplot(aes(x = GRMZM2G046776, y = GRMZM2G125969, color = Type)) + 
  geom_point(size = I(3), alpha = 1/2) + 
  scale_colour_brewer(type = "qual", palette = "Dark2") + 
  theme(legend.position = "none")

sca_rna2 <- dd_sc  %>%
  dplyr::select(imp_genes10,Type ) %>%
  ggplot(aes(x = GRMZM2G093708, y = GRMZM2G014695, color = Type)) +
  geom_point(size = I(3), alpha = 1/2) + 
  scale_colour_brewer(type = "qual", palette = "Dark2") + 
  theme(legend.position = "none")

sca_rna3 <- dd_sc  %>%
  dplyr::select(imp_genes10,Type ) %>%
  ggplot(aes(x = GRMZM2G125969, y = GRMZM2G075892, color=Type)) +
  geom_point(size = I(3), alpha = 1/2) + 
  scale_colour_brewer(type = "qual", palette = "Dark2") +
  theme(legend.position="none")

ppfimpo1_rna <- ggplot( globalimpo_rna[1:50,], aes(x = mean, y = variable)) + 
  geom_point() + 
  theme(aspect.ratio = 1, axis.text.y = element_blank()) +
  labs(x = "Global average importance", y="")

(pl_rna1|ppfimpo1_rna)/(sca_rna1|sca_rna2|sca_rna3 )

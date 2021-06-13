#install.packages('shiny')
#install.packages('shinythemes)
#install.packages('mixtools')
library(shiny)
library(shinythemes)
library(mixtools)

ui <- fluidPage(
    theme=shinytheme("slate"),
    titlePanel("Approximation de toute densité continue par un mélange gaussien"),
    mainPanel(width = 12,
        h2("Modèle de mélange gaussien"),
        p('Un modèle de mélange gaussien est un modèle statistique exprimé par une densité de mélange. Il sert usuellement à estimer paramétriquement la distribution de variables aléatoires en les modélisant comme une somme de plusieurs gaussiennes.'),
        h2("Le modèle"),
        p('Le modèle de mélange fini de lois de probabilité consiste à supposer que les données proviennent d’une source contenant plusieurs sous-populations homogènes appelées composants.
        La population totale est un mélange de ces sous-populations. Le modèle résultant est un modèle de mélange fini.'),
        withMathJax(p("Soit \\(X=(X_1,...,X_n)\\) un échantillon de variables aléatoires indépendantes identiquement distribuées (iid ) de loi de mélange fini à \\(K\\) composants, de densité f dont  la forme générale est :
        $$f(x)=\\sum_{k=1}^{K}\\pi_k h_{\\Phi_k}(x))$$ avec:"),
        withMathJax(p("-\\(\\pi_k\\) 
        les proportions respectives des sous populations telle que \\(0<\\pi_k\\)\\(\\leq1\\) et \\(\\sum_{k=1}^{K}\\pi_k=1\\) ")),
        withMathJax(p("-\\(h_{\\Phi_k}\\) la densité du k-ième composant (une gaussienne \\(\\mathcal{N}(\\mu_k,\\sigma_k^2)\\))")),
        withMathJax(p("Le modèle de mélange est un modèle à données manquantes. En effet, si on échantillonnait dans une population formée par K sous-populations, on devrait avoir les couples (\\(X_i,Z_i\\)), où \\(X_i=x_i\\) indique la mesure faite sur le i-ème individu et \\(Z_i=k\\) indique 
        le numéro de la sous-population à laquelle appartient cet individu. 
        En échantillonnant dans la sous population \\(k\\) et en supposant \\(X\\) discrète, on obtiendrait alors le modèle \\(\\mathbb{P}( X =x|Z=k )= f_k( x ,\\alpha_k)\\). Mais le paramètre \\(\\alpha_k\\) est en général inconnu et propre à la k-ième sous-population.
        De même,les données manquantes sont \\(Z=(Z_1,\\dots, Z_n)\\),avec \\(Z_i=k\\) si \\(i\\) provient du groupe \\(k\\).
        On n'observe donc que l’échantillon \\(( X_1,\\dots,X_n)\\). 
        Le mélange définit plus haut peut être vu comme la loi marginale de la variable \\(X\\) pour le couple \\((X , Z)\\). C'est donc un modèle à données manquantes.
        ")),
        h3("Remarque"),
        h4("Toute densité continue peut s'approcher à l'aide d'un mélange gaussien, au sens de la norme \\(L_1\\) ou uniformement sur tout compact."),
        p("Cela s'appuie sur le théorème suivant:"),
        h2("Théorème"),
        withMathJax(p("Soit \\(g\\) une densité continue. Pour tout \\(\\varepsilon >0\\) il existe un mélange gaussien fini de densité \\(\\overline{g}\\)  donnée par :
        $$\\overline{g}=\\sum_{j=1}^{m} \\pi_jf_\\mathcal{N(\\mu_j,\\sigma_j^2)}(x) $$
        avec \\(\\mu_j\\in \\mathbb{R}\\), \\(\\sigma_j >0\\) tel que \\(\\sum_{j=1}^{m}\\pi_j=1\\), tel que 
        $$\\|g-\\overline{g}\\|_1:=\\int_\\mathbb{R} \\ |\\ g(x)-\\overline{g}(x)|<\\varepsilon $$
        De même, pour tout \\(\\varepsilon >0\\) et tout compact \\(\\kappa\\), il existe un mélange gaussien fini de densité \\(\\overline{g}\\) tel que
        $$\\sup_{x\\in \\kappa} |\\ g(x)-\\overline{g}(x) |<\\varepsilon $$")),
        h2("Illustration du théorème"),
        withMathJax(p("Ici, notre but sera d'illustrer graphiquement ce thèorème. On se restreindra aux densités continues usuelles: densité d'une loi gamma de paramètres \\(a\\) et \\(b\\) , \\(\\Gamma(a,b)\\), loi Beta \\(B(a,b)\\), loi de Student \\(T(k)\\) à \\(k\\) degrés de liberté, loi du \\(\\mathcal{\\chi}^2_{(d1)}\\) et loi de Fisher \\(F(d_1,d_2)\\). 
        On verra par la suite que toutes ces densités peuvent être approchées par des mélanges gaussiens à densité unimodale, parfois dyssimétrique ou fortement dyssimétrique, parfois par des densités avec outlier. Comment approcher donc ces densités?
        En reprenant la formule plus haut, $$\\overline{g}=\\sum_{j=1}^{m} \\pi_jf_\\mathcal{N(\\mu_j,\\sigma_j^2)}(x) $$
        on remarque que pour pouvoir déterminer notre mélange gaussien, on devra estimer le vecteur des paramètres \\(\\theta = (\\mu_1,...,\\mu_m, \\sigma_1^2,...,\\sigma_m^2,\\pi_1,...,\\pi_m)\\), \\(j \\in \\ [1,..m]\\).
        Pour le faire, on a décidé d'utiliser l'algorithme EM (espérance-maximisation). Cet algorithme calcule par itération les paramètres de l'EMV en considérant l'esperance conditionnelle de la log-vraisemblance des données observées et le paramètre obtenu à l'itération précedente.
        On a décidé d'utiliser cet algorithme à travers le package R 'mixtools' (pour les 3 premiers exemples ) et l'algorithme EM vu en TP (pour les 2 derniers exemples ).")),
        h2("Simulations"),
        withMathJax(p("Passons maintenant à la simulation pour les différentes lois! Vous pouvez choisir les paramètres suivants:")),
        withMathJax(p("- \\(n\\): nombre d'observations qui vont être utilisées par l'algorithme")),
        withMathJax(p("- \\(d_1,  d_2 , a, b\\): paramètres des lois choisies")),
        withMathJax(p("- \\(l\\): ordre du mélange ")),
        sidebarPanel(
        selectInput("model",label="Loi:",choices=c("Chi 2","Student","Fisher","Beta","Gamma"),selected="Chi 2"),
        numericInput("n",label="Nombre d'observations:",min = 100,max = 10000,value = 1000),
        conditionalPanel(condition = "input.model == 'Chi 2'",
                         numericInput("l",label="Ordre du mélange:",min=2,max=4,value = 2),
                         sliderInput("k",label="Degrés de liberté :",min = 0,max = 20,value = 8, step=1)),
        conditionalPanel(condition = "input.model == 'Student'",
                         numericInput("l",label="Ordre du mélange:",min=2,max=4,value = 2),
                         sliderInput("d",label="Degrés de liberté :",min = 0,max = 20,value = 8, step=1)),
        conditionalPanel(condition="input.model=='Fisher'",
                         numericInput("l",label="Ordre du mélange:",min=2,max=4,value = 2),                   
                         sliderInput("df1",label="Degrés de liberté 1:",min = 0,max = 20,value = 8, step=1),
                         sliderInput("df2",label="Degrés de liberté 2:",min = 0,max = 10,value = 6, step=1)),
        conditionalPanel(condition="input.model=='Beta'",
                         sliderInput("alpha",label="Paramètre alpha",min=0.5,max=5,value=2,step=0.5),
                         sliderInput("beta",label="Paramètre beta",min=0.5,max=5,value=3,step=0.5)),
        conditionalPanel(condition="input.model=='Gamma'",
                         sliderInput("a",label="Paramètre a",min=0.5,max=5,value=2,step=0.5),
                         sliderInput("b",label="Paramètre b",min=0.5,max=5,value=2,step=0.5))),
        mainPanel(plotOutput("studPlot"))
        
        )
            
)

)

    


server <- function(input, output) {
    
    dchi2 <- function(x,k){
        dens <- rep(0,length(x))
        return(dens + dchisq(x,k))
        
    }
    
    
    dstudent <- function(x,d){
        dens <- rep(0, length(x))
        return (dens+dt(x,d))
    }
    
    
    dfisher <- function(x,df1,df2){
        dens <- rep(0,length(x))
        return (dens +df(x,df1,df2))
    }
    
    
    melange <- function(x,mix){
        dens <- rep(0, length(x))
        K <- length(mix$lambda)
        for (k in 1:K){
            dens <- dens + mix$lambda[k]*dnorm(x,mix$mu[k],mix$sigma[k])
            return (dens)
        }
    }
    
    
    update.theta <- function(obs, theta){
        nb <- length(obs) 
        K  <- length(theta$mu) 
        
        # calcul des alpha_ji
        alpha <- matrix(NA, nrow=K, ncol=nb)
        for (j in 1:K){
            alpha[j,] <- theta$pi[j]*dnorm(obs,theta$mu[j],theta$sig[j])
        }
        p.theta <- apply(alpha, 2, 'sum')
        alpha   <- alpha/matrix(p.theta, byrow=TRUE, nrow=K, ncol=nb)
        pi.new  <- apply(alpha, 1, 'mean')
        mu.new  <- c(alpha%*%obs/nb/pi.new)
        mat     <- (matrix(rep(obs,each=K), nrow=K) -matrix(rep(mu.new, times=nb), nrow=K))^2
        sig.new <- sqrt(apply(alpha*mat, 1, 'mean')/pi.new)
        
        theta.new <-  list(pi =pi.new, mu=mu.new, sig=sig.new)
        return(theta.new)
    }
    
    
    rnormmix <- function(n, theta){
        K <- length(theta$pi)
        etiqu <- sample(1:K, size=n, replace=TRUE, prob=theta$pi)
        obs <- rep(NA, n)
        i <- 1
        for (k in 1:K){
            n.group <- sum(etiqu==k)
            if(n.group>0){
                obs[i:(i+n.group-1)] <- rnorm(n.group, theta$mu[k], theta$sig[k])
                i <- i + n.group
            }
        }        
        return(obs)
    }
    
    
    algoEM <- function(obs, theta.init){   
        
        resultat <- list(emv = theta.new, nb.iter = it)
        return(resultat)
    }
    
    
    dnormmix <- function(x,theta){
        dens <- rep(0, length(x))
        K <- length(theta$pi)
        for (k in 1:K){
            dens <- dens + theta$pi[k]*dnorm(x, theta$mu[k], theta$sig[k])
        }                                         
        return(dens)
    }
    
    lvraisnorm <- function(param, obs){    
        logvrais <- sum(log(dnormmix(obs, param)))
        return(logvrais)
    }
    
    
    algoEM <- function(obs, theta.init, R=200, epsilon=1e-3){   
        theta.old <- theta.init
        crit_arret <- FALSE
        log.vrais.old <- lvraisnorm(theta.init, obs)
        it <- 0
        while (!crit_arret && (it < R))
        { 
            theta.new <- update.theta(obs, theta.old)
            log.vrais.new <- lvraisnorm(theta.new, obs)
            crit_arret <- (abs((log.vrais.new - log.vrais.old)/log.vrais.old) < epsilon)
            log.vrais.old <- log.vrais.new
            theta.old <- theta.new
            it <- it + 1
        }
        resultat <- list(emv = theta.new, nb.iter = it)
        return(resultat)
    }
    
    
    
    output$studPlot <- renderPlot({
        
        n <- input$n
        l <- input$l
        k <- input$k
        d <- input$d
        d1 <- input$df1
        d2 <- input$df2
        a<-input$alpha
        b<-input$beta
        a1<-input$a
        b1<-input$b
        
        
        if (input$model=='Chi 2'){
            
            obs <- rchisq(n,k)
            hist(obs,xlim=c(0,20),ylim=c(0,.4),prob =TRUE,main=NULL,lty=3)
            curve(dchi2(x,k),xlim=c(0,20),ylim=c(0,.4),lwd=4,ylab='densité',xlab='observations',col='pink',add = TRUE)
            
            mix <- normalmixEM(obs,k=l,maxit=100,epsilon=.01)
            
            
            curve(melange(x,mix),xlim=c(0,20),ylim=c(0,.4),add=TRUE,col='red',lwd=4)
            legend(15,.3,c('Densité Chi-2','Mélange gaussien'),col=c('pink','red'),pch=16)
            
        }
        
        
        if (input$model=="Student"){
            
            obs <- rt(n,d)
            hist(obs,xlim=c(-4,4),ylim=c(0,.4),prob =TRUE,main=NULL,lty=3)
            curve(dstudent(x,d),xlim=c(-4,4),ylim=c(0,.6),lwd=4,ylab='Densité',xlab='observations',col='green',add = TRUE)
            
            mix <- normalmixEM(obs,k=l,maxit=100,epsilon=.01)
            
            
            curve(melange(x,mix),xlim=c(-4,4),ylim=c(0,.5),add=TRUE,col='red',lwd=4)
            legend(2.5,.3,c('Densité Student','Mélange gaussien'),col=c('green','red'),pch=16)
        }
        
        
        if (input$model=='Fisher'){
            
            obs <- rf(n,d1,d2)
            hist(obs,xlim=c(0,11),ylim=c(0,.7),prob =TRUE,main=NULL,lty=3)
            curve(dfisher(x,d1,d2),xlim=c(0,11),ylim=c(0,.7),lwd=4,ylab='Densité',xlab='Observations',col='orange',add = TRUE)
            
            mix <- normalmixEM(obs,k=l,maxit=100,epsilon=.01)
            
            
            curve(melange(x,mix),xlim=c(0,11),ylim=c(0,.7),add=TRUE,col='red',lwd=4)
            legend(7,.3,c('Densité Fisher','Mélange gaussien'),col=c('orange','red'),pch=16)
        }
        
        
        if (input$model=='Beta'){
            
            obs <- rbeta(n,a,b)
            theta = list(pi=c(1,1,3)/5, mu=c(0,1/2,13/15), sig=c(1,(2/3)^2,(5/9)^2))
            hist(obs,xlim=c(0,1.2),ylim=c(0,2.5),prob =TRUE,main=NULL,lty=3)
            curve(dbeta(x,a,b),lwd=4,ylab='densité',xlab='observations',col='blue',add = TRUE,xlim=c(0,1.2),ylim=c(0,2.5))
            
            res <- algoEM(obs,theta)
            
            curve(dnormmix(x, res$emv), add=TRUE, col='red')
            legend(1,2,c('Densité Beta','Mélange gaussien'),col=c('blue','red'),pch=16)
            
        }
        
        
        if (input$model=='Gamma'){
            
            obs2 <- rgamma(n,a1,b1)
            theta2 = list(pi=c(1,1,3)/5, mu=c(0,1/2,13/15), sig=c(1,(2/3)^2,(5/9)^2))  
            hist(obs2,xlim=c(0,6),ylim=c(0,1.2),prob =TRUE,main=NULL,lty=3)
            curve(dgamma(x,a1,b1),lwd=4,ylab='densité',xlab='observations',col='blue',add = TRUE)
            
            res2 <- algoEM(obs2,theta2)
            
            curve(dnormmix(x, res2$emv), add=TRUE, col='red')
            legend(4,1,c('Densité gamma','Mélange gaussien'),col=c('blue','red'),pch=16)
            
        }   
        
        
    })
    
}
# Run the application 
shinyApp(ui = ui, server = server)
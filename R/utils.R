# Series of util functions to evaluate the MCMC chain for the bpwpm model

plot_run <- function(mcmc_chain, plot_dist = c(2,2), n = 10^3){
    " Function that plots each one of the chains on a separate plot.

    INPUT:
        mcmc_chain:= Matrix of n*p, where n is the number of draws and p is the number of parameters
        plot_dist:= Distribution of the plots on the graphical interphase

    "
    # Graphical parameters
    par(mfrow = plot_dist)

    dim_mcmc <- dim(mcmc_chain)
    draws <- dim_mcmc[1]
    params <- dim_mcmc[2]

    n <- min(n,draws)

    for(i in 1:params){

        plot(x = 1:n, y = mcmc_chain[1:n, i], type = "l")
        lines(mcmc_chain[1:n,i])

        # if((i %% prod(plot_dist)) == 0){
        #     dev.off()
        #     readline(prompt = "Pause. Press <Enter> to continue...")
        # }
    }
}

#-------------------------------------------------------------------------------

thin_chain <- function(x){
    return(x)
}

#-------------------------------------------------------------------------------

# library(RColorBrewer)
# analisis_rapido <- function(modelo, draws = 1000,
#                             summary = FALSE, each = FALSE){
#     # Función que analiza la convergencia de forma muy general de un modelo
#
#     # Resumen general del modelo
#     if(summary == TRUE){
#         print(summary(modelo))
#     }
#
#     # Variables
#     nombre_vars <- dimnames(modelo$param)[[2]]
#     n <- dim(modelo$param)
#
#     # Tomamos la última parte de la cadena para que no se vea tanto
#     vars_MCMC <- modelo$param[(n[1]-draws):n[1], ]
#     colnames(vars_MCMC) <- nombre_vars
#
#     # Pendejada para sacar colores y que la gráfica tenga sentido
#     qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#     col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
#                                rownames(qual_col_pals)))
#
#     # Matriz con todas las cadenas
#     matplot(1:dim(vars_MCMC)[1], vars_MCMC, type = 'l',
#             xlab  = "Cadena", ylab = "Beta", col = col_vector)
#     legend('right', legend = colnames(vars_MCMC), cex = .5, pch = 0,
#            fill = col_vector)
#
#     # Por si quieres ver cada una
#     if(each == TRUE){
#         for(i in 1:dim(vars_MCMC)[2]){
#             variable <- i # Cambiar este número para ver las diferentes variables
#             plot(modelo_coahuila$param[,variable], type = 'l',
#                  main = nombre_vars[variable], col = "red")
#             readline(prompt = "Pause. Press <Enter> to continue...")
#         }
#     }
# }

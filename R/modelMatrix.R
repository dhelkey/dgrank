                                             
modelMatrix = function(x, interactions = FALSE, sparse = TRUE, intercept = FALSE){
    #Constructs a model.matrix
    #Factors use treatment encoding, i.e. indicators for 2nd through nth level
    #TODO Currently has problems with just 1 continious variable losing name information
    #TODO Currently has problems with just 1 continious variable losing name information
    options(contrasts = rep("contr.treatment", 2))
    model_formula = as.formula( ~.) #First order Effects
    if (interactions){model_formula = as.formula( ~(.)^2)} #1st order and interctions
    mat = model.matrix(model_formula, data = data.frame(x))
    if (!intercept){mat = mat[ ,-1]}
    return( as.matrix(mat) )
}
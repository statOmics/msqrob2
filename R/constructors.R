msqrobModel<-function(modelType,model)
{
mod<-new("msqrobModel")
mod@modelType<-modelType
mod@model<-model
return(mod)
}

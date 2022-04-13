list.initialize = function (var.len, ini.value=NULL) {

ans=list()
for (m in 1:var.len) {
    ans[[m]]=ini.value
}
return(ans)

}
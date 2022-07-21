# socialmixr equation
#  C'ij = 1/2Ni[CijNi+CjiNj]
#

rc <- matrix(c(20,10,5,10,15,10,5,10,5),nrow = 3)
rownames(rc) <- c("Y","A","E")
colnames(rc) <- c("Y","A","E")

N  <- c(Y=3000,A=5000,E=2000)
Tr <- c(Y=0.02,A=0.04,E=0.04)
h  <- 0.5

c <- rc
c[2,1] <- 1/(2*N[2])*(rc[2,1]*N[2]+rc[1,2]*N[1])
c[1,2] <- 1/(2*N[1])*(rc[1,2]*N[1]+rc[2,1]*N[2])
c[3,1] <- 1/(2*N[3])*(rc[3,1]*N[3]+rc[1,3]*N[1])
c[1,3] <- 1/(2*N[1])*(rc[1,3]*N[1]+rc[3,1]*N[3])
c[2,3] <- 1/(2*N[2])*(rc[2,3]*N[2]+rc[3,2]*N[3])
c[3,2] <- 1/(2*N[3])*(rc[3,2]*N[3]+rc[2,3]*N[2])

tc <- N*c

IS <- c(Y=0,A=0,E=0)
IP <- c(Y=0,A=0,E=0)
IC <- c(Y=0,A=10,E=0)

# FOI detailed
lambda1 <- matrix(rep(0,9),nrow = 3)
rownames(lambda1) <- c("Y","A", "E")
colnames(lambda1) <- c("Y","A", "E")

lambda1["Y","Y"] <- Tr["Y"] * c["Y","Y"] * (h*IS["Y"] + IP["Y"]+IC["Y"])/N["Y"]
lambda1["Y","A"] <- Tr["Y"] * c["Y","A"] * (h*IS["A"] + IP["A"]+IC["A"])/N["A"]
lambda1["Y","E"] <- Tr["Y"] * c["Y","E"] * (h*IS["E"] + IP["E"]+IC["E"])/N["E"]

lambda1["A","Y"] <- Tr["A"] * c["A","Y"] * (h*IS["Y"] + IP["Y"]+IC["Y"])/N["Y"]
lambda1["A","A"] <- Tr["A"] * c["A","A"] * (h*IS["A"] + IP["A"]+IC["A"])/N["A"]
lambda1["A","E"] <- Tr["A"] * c["A","E"] * (h*IS["E"] + IP["E"]+IC["E"])/N["E"]

lambda1["E","Y"] <- Tr["E"] * c["E","Y"] * (h*IS["Y"] + IP["Y"]+IC["Y"])/N["Y"]
lambda1["E","A"] <- Tr["E"] * c["E","A"] * (h*IS["A"] + IP["A"]+IC["A"])/N["A"]
lambda1["E","E"] <- Tr["E"] * c["E","E"] * (h*IS["E"] + IP["E"]+IC["E"])/N["E"]

foi_lambda1 <- apply(lambda1,1,sum)

lambda2 <- matrix(rep(0,9),nrow = 3)
rownames(lambda2) <- c("Y","A", "E")
colnames(lambda2) <- c("Y","A", "E")

for(i in seq_along(c("Y","A","E"))){
  for(j in seq_along(c("Y","A","E"))){
    lambda2[i,j] <- Tr[i] * c[i,j] * (h*IS[j]+IP[j]+IC[j])/N[j]
  }
}

foi_lambda2 <- apply(lambda2,1,sum)

test1 <- lambda1 - lambda2
test2 <- foi_lambda1 - foi_lambda2






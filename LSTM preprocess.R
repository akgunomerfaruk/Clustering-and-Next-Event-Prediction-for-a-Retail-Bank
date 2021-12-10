library(dplyr)
library(tidyr)
library(caret)
library(tensorflow)
library(keras)

set.seed(499)

freqnew <- read.csv(file = 'sequencedata.csv', header=FALSE)
freqnew <-freqnew[apply(freqnew[,-1], 1, function(x) !all(x==0)),]
bankevents <- read.csv(file = 'bank events.csv', header=TRUE)
bankevents <- as.character(bankevents[ ,1])

freqnew <- data.matrix(freqnew[,-1])
seq <- Map(Filter, list(Negate(is.na)), split(freqnew, 1:nrow(freqnew)))

padseq <- pad_sequences(seq, maxlen = 3456, dtype = "int32",
              padding = "pre", truncating = "pre", value = 0)

customerseq <- lapply(seq, function(A) A[A %in% bankevents == FALSE])

padcust <- pad_sequences(customerseq, maxlen = 1579, dtype = "int32",
                         padding = "pre", truncating = "pre", value = 0)

padcust <-padcust[apply(padcust[,-1], 1, function(x) !all(x==0)),]

X1 <- padcust[ ,1:1578]
Y1 <- as.factor(padcust[ ,1579])
Y <- data.frame(id=c(1:length(Y1)), event=Y1)

dmy <- dummyVars(" ~ .", data = Y)
Y_hot <- data.frame(predict(dmy, newdata = Y))
Y_hot <- Y_hot[,-1]

reshape_X_3d <- function(X) {
  dim(X) <- c(dim(X)[1], dim(X)[2], 1)
  X
}

X <- reshape_X_3d(X1)

write.csv(X1, 'X.csv')
write.csv(Y_hot, 'Y.csv')

a <- read.csv(file = 'X.csv', header=TRUE)
b<-  read.csv(file = 'Y.csv', header=TRUE)

model <- keras_model_sequential() %>%   
  layer_lstm(units=256, input_shape = c(1578,1)) %>%   
  layer_dropout(0.2) %>%
  layer_dense(units=128, activation= "relu") %>%
  layer_dropout(0.2) %>%
  layer_dense(units=58, activation = "softmax")

model %>% compile(loss = 'categorical_crossentropy',
                  optimizer = 'adam',
                  metrics = 'auc'
)

remove.packages("tensorflow")
model %>% summary()

model %>% fit(X,Y_hot, epochs=75, batch_size=10, validation_split=0.2, shuffle=FALSE)
y_pred <- model %>% predict(X)

scores = model %>% evaluate(X, y, verbose = 0)
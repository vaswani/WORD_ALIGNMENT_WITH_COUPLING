#! /bin/bash

../GIZA++-v2/GIZA++ -c train.en_train.fr.snt -coocurrencefile train.en_train.fr.cooc -s train.en.vcb -t train.fr.vcb -cef train.fr_train.en.snt -efcoocurrencefile train.fr_train.en.cooc -cfe train.en_train.fr.snt -fecoocurrencefile train.en_train.fr.cooc -model1iterations 4 -model2iterations 0 -model3iterations 0 -model4iterations 0 -model5iterations 0 -model6iterations 0 -mh 0 -model1dumpfrequency 1 -hmmdumpfrequency 1 -eta 0.9 -regLambda 10000 --joint_reg 1 --joint_reg 0

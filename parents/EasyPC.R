 
###################################################################################
#                              Easy Pedigree Checking (EasyPC)                    #
# A program for paternity test(biparental && single ) with high density SNP data. #
#                                         Designed by Zhe Zhang,Yuanyu Luo et al. #    
#                                                     Email: zhezhang@scau.edu.cn #
###################################################################################

#####################################################################################
#                   A brief introduction for the program
# 该程序主要利用高密度的SNP芯片数据进行亲子鉴定，双亲和单亲的亲子鉴定均可。调用方�?: ParTest(arg1, arg2,...)
# 输入参数释义�?
#           'ped'为系谱数�?
#           'gen'表示个体的基因型数据�?
#           'thre'代表阈�?
#           'is.picture'设置是否生成孟德尔错误率的直方图 
#           'is.merge.output'设置write到文件夹的错误系谱是否包含推断出的正确ID
#           'stop.if.snp.null'表示是否系谱中出现没有基因型的个体就终止程序   
#  输出结果�?
#         (1)一个命名为'paternity_test(Biparental or Single)'文件夹，包含以下�?
#           'ped_check_result':文件夹，主要调用PedCheck()函数得到的错误和正确的系谱记录。一般包�?'ped_right_record.txt'�?'ped_wrong_record.txt'等文�?
#           'ped_infer_result':文件夹，主要调用Pedinfer()函数得到的推断出的正确系谱记录。一般包�?'ped_right_infer.txt', 'sirID_wrong_record.txt',
#                              'damID_wrong_record.txt'�? 'offID_wrong_record.txt'等文件。注意如果没有就不会生成，例如系谱中没有sirID错误的系谱记录就不会生成'sirID_wrong_record.txt'
#           'err_bargraph_01':png格式图，孟德尔错误率的直方图，包含孟德尔错误率小于阈值的情况，及所有系谱记录的错误率直方图
#           'err_bargraph_02':png格式图，孟德尔错误率的直方图, 但不包含孟德尔错误率小于阈值的情况，及只有错误系谱记录的错误率直方�?
#           'mend_err_res':'txt'文件，所有系谱记录的错误数和错误�?
#         (2)R console上会显示有相关结果和提示内容，主要包含以下：
#           'Res_01:SNP Date Check':调用SnpDataCheck()函数检查系谱中是否有缺失基因型的个体，并进行相应提�?
#           'Res_02:Ped Check':主要调用PedCheck()函数对系谱中是否有错误系谱记录进行检查，并进行相应提�?
#           'Res_03:Ped Infer':主要调用PedInfer()函数帮助存在于错误系谱记录中的子代个体重新找回其亲本，并进行相应输出
# 注意事项�?
#         系谱命名必须规范：对于单亲，第一列名�?'offID'，第二列'parID';对于双亲，第一列名�?'offID', 第二�?'sirID', 第三�?'damID'
#         'gen'表示符规范：'gen'第一列是ID号；ID号后面的内容为SNP位点基因型信息，只能�?'1', '2', '3', '5'表示，其�?'1'表示"aa", '2'表示'Aa', '3'表示'AA', '5'代表'unknown'
#################################################################################################################################################################################################



ParTest <- function(ped, gen, thre=0.01, is.picture=T, is.merge.output=F, stop.if.snp.null=T, ...){
    # The following function is designed for checking the individuals without SNP data
    SnpDataCheck <- function(ped, gen, ...){
        ped.colnames <- colnames(ped)   
        cat("Res_01:SNP Data Check\n") 
        cat("(1)Individuals' SNP data checked as follows:\n") 
        tmp1 <- c()  # 存储没有SNP数据的ID号下标的向量
        for(i in 1:ncol(ped)){    
            t <- 0
            ped.only.id <- unique(ped[, i])  # 系谱中由于存在一个ID号同时出现多次，所以取系谱中ID号的唯一值，便于检�?
            if(any(tmp <- !ped.only.id %in% gen[, 1])){  # 系谱中ID号没有出现在gen中，说明SNP数据缺失，注意这里只是记录TRUE or FALSE的数据类�?
                t <- 1
                id.without.snp <- ped.only.id[tmp]  # 缺失SNP数据的ID�?
                if(!stop.if.snp.null){  # 系谱中个体SNP数据缺失�? 如果程序需要继续执行，则要将包含这些ID号的记录剔除
                    tmp2 <- match(id.without.snp, ped[, i])  # 记录缺失SNP的个体ID号分别在系谱中的位置
                    tmp1 <- c(tmp1, tmp2)  # 记录所有缺失SNP数据个体ID号的位置信息
                }
                cat(paste0("[", i, "]"), "The", ped.colnames[i], "without the SNP data:", paste(sQuote(id.without.snp), collapse=","), "\n     the total:", length(id.without.snp), sep=" ", "\n") 
            }else{
                cat(paste0("[", i, "]"), "All the", ped.colnames[i], "have the SNP data", fill=TRUE, sep=" ","\n") 
            }     
        } 
        if(t){  # 如果t==0�? 则说明系谱中不存在缺失SNP的个体ID，至此该函数可以结束。否则，按照设置的参数要求继续执�?
            if(stop.if.snp.null){  # 如果stop.if.snp.null==T, 则程序结束，并输出提�?
                stop("There are some individuals without snp data, please check it!", call.=FALSE)
            }else{  # 否则�? 将包含缺失SNP数据的ID号系谱记录剔除后，程序继续执�?
                record.del <- ped[nr <- unique(tmp1), ]
                cat("(2)The pedigree", ngettext(length(nr), "record is", "records are"), "deleted for snp data missing\n")
                print(record.del)
                ped <- ped[-nr, ]
            }
        }
        return(ped)
    }

    # The following function is designed for typecasting, 'data.frame' ———�?>> 'matrix'
    TypeCast <- function(x){  
        x <- matrix(as.numeric(as.matrix(x)), ncol=ncol(x))  # 类型转换函数�? 将数据框转化为矩阵形�?
        return(x)
    }

    # The following function is designed for Mendelian Error&&Errate  
    MendErrCal <- function(off.gen, par.one.gen, par.two.gen=NULL){  # 参数par.two.gen默认为NULL，表示为单亲的情�?
        if(is.null(par.two.gen)){  # 单亲情况下的孟德尔错误和错误率的计算方式
            if(any(c(off.gen, par.one.gen) == 5)){ # '5'表示SNP位点unknown的情况， 如果位点缺失则视为无效位�?
                tmp <- which((off.gen*par.one.gen) %% 5 == 0) 
                off.gen <- off.gen[-tmp]
                par.one.gen <- par.one.gen[-tmp]
            }
            gen.dev <- off.gen/par.one.gen  # SNP位点比，用于错误位点的筛�?
            num.of.err <- sum(c(gen.dev == 3) | c(gen.dev == 1/3))  # 如果亲本�?'3':AA, 子代'1':aa，则不符合孟德尔遗传定律；反之亦然�?
            errate <- num.of.err / sum(par.one.gen != 2)  # 计算错误率时排除亲本基因�?'2':Aa的情�?
            return(list(num.of.err=num.of.err, errate=errate))
        }else{
            if(any(c(off.gen, par.one.gen, par.two.gen) == 5)){  # 同上单亲情况分析   
                    tmp1 <- which((off.gen*par.one.gen*par.two.gen) %% 5 == 0) 
                    off.gen <- off.gen[-tmp1]
                    par.one.gen <- par.one.gen[-tmp1]
                    par.two.gen <- par.two.gen[-tmp1]  
            }
            tmp2 <- which(par.one.gen*par.two.gen == 4)  # 双亲SNP位点�?'2':Aa的情况无法判�?
            off.gen <- off.gen[-tmp2]
            par.one.gen <- par.one.gen[-tmp2]
            par.two.gen <- par.two.gen[-tmp2]
            gen.dev1 <- off.gen/par.one.gen  # 同上单亲情况分析
            gen.dev2 <- off.gen/par.two.gen
            gen.dev3 <- gen.dev1/par.two.gen
            num.of.err <- sum(c(gen.dev1 == 3) | c(gen.dev2 == 3) | c(gen.dev3 == 2) | c(gen.dev1 == 1/3) | c(gen.dev2 == 1/3) | c(gen.dev3 == 2/9))- sum(c(gen.dev3 == 3) | c(gen.dev3 == 1/9))  # 根据孟德尔定律，双亲情况下共有六种错误情况，两种重复排除情况
            errate <- num.of.err/length(off.gen) 
            return(list(num.of.err = num.of.err, errate = errate))      
        }
    }

    # The following function is designed for drawing the histogram of errate
    ErratePlot <- function(errate, nc){  # 参数nc表示SNP位点�?
        hist(as.numeric(errate), 
             main = "Paternity Test\n\n",
             xlab = "Errate", 
             ylab = "Frequency", 
             breaks = length(errate)/2       
             )
        title(main=paste0("(SNP=", nc, ")"), sub=date(), cex.sub=.5, cex.main=.95)
    }

    # The following function is designed for searching for the right parentage in the pedigree
    ParentInfer <- function(ped, gen, ped.wrong.record, ped.col=2, ...){
        if(nrow(gen) > length(ped.inds <- unique(as.vector(as.matrix(ped))))){  # 如果提供gen的个体比系谱中个体数多， 则推断时需要以gen中个体数为基础
            par.id.set <- c(unique(ped[, ped.col]), sort(gen[!gen[, 1] %in% ped.inds, 1], decreasing = as.logical(ped.col-2)))  # 系谱中个体和gen集中不包含在ped中的个体构成亲本ID集合
        }else{
            par.id.set <- unique(ped[, ped.col])
        }
        tmp1 <- tmp2 <- matrix(, 0, 2, dimnames=list(NULL, c("offID", ngettext(ped.col-1, "sirID", "damID"))))   
        par.id.index <- seq_along(par.id.set)  # 记录亲本ID的下�?
        off.id.set <- ped.wrong.record[, 1]  #  错误系谱记录的子代ID集合
        off.gen.set <- TypeCast(gen[match(off.id.set, gen[, 1]), -1])  # 与子代ID号对应的子代SNP位点数据�?
        par.gen.set <- TypeCast(gen[match(par.id.set, gen[, 1]), -1])  # 同上分析
        for (i in seq_along(off.id.set)) {
            if(mode == 2){  # 'mode==2'表示双亲情况
                id.per.seq <- c(which(par.id.set == ped.wrong.record[i, ped.col]), par.id.index[!par.id.set %in% ped.wrong.record[i, 1]])  # 为提高亲本推断的效率�? 亲本ID号按照一定规则重�?      
                row.index <- 0
            }else{
                id.per.seq <- par.id.index[-which(par.id.set == ped.wrong.record[i, 2])] 
            }
            par.gen.tmp <- par.gen.set[id.per.seq, ]  # 亲本ID重排后对应的亲本gen�?
            par.off.pair <- expand.grid(off.id = off.id.set[i], par.id.set = par.id.set[id.per.seq]) 
            for (j in seq_along(id.per.seq)) {
                if(ncol(off.gen.set) <= 10000){  # 设置if条件语句是为提高推断效率采取的一种策�?
                    if(MendErrCal(off.gen.set[i, ], par.gen.tmp[j, ])[[2]] <= thre){
                        row.index <- j
                        break       
                    }                    
                }else{
                    if(MendErrCal(off.gen.set[i, 1:1000], par.gen.tmp[j, 1:1000])[[2]] <= thre) {
                        if(MendErrCal(off.gen.set[i, ], par.gen.tmp[j, ])[[2]] <= thre){
                            row.index <- j
                            break       
                        }
                    }                
                }  
            }
            if(mode == 2 & row.index == 1){
                tmp1 <- rbind(tmp1, as.matrix(par.off.pair[row.index, ]))  # 双亲情况下如�?'row.index==1'�? 则说明原系谱中的"off_sir",�?"off_dam"是正确的�? 因为是进行单亲检测所以排除性的断定是另一亲本记录�?
            }
            tmp2 <- rbind(tmp2, as.matrix(par.off.pair[row.index, ]))  # 推断出的正确亲子关系，注意一次只能推断出一个亲�?
        }
        if(mode == 2){
            tmp3 <- ped.wrong.record[match(tmp1[, 1], ped.wrong.record[, 1]), ]   
            tmp3 <- tmp3[order(as.numeric(tmp3[, 1])), ]
            rownames(tmp3) <- NULL
            return(list(tmp2=tmp2, tmp3=tmp3))               
        }else{
            colnames(tmp2) <- c("offID", "parID")
            return(tmp2)
        }
    }

    # The following function is designed for paternity(single) test 
    SinPedCheck <- function(ped, gen, ...){
        dir.create("ped_check_result")  # 结果输出�?'ped_check_result'文件夹中
        setwd("./ped_check_result")
        off.id <- ped[, 1]
        par.id <- ped[, 2]
        off.gen.set <- gen[match(off.id, gen[, 1]), -1]
        par.gen.set <- gen[match(par.id, gen[, 1]), -1]
        nr<-nrow(ped)
        mend.err.res <- data.frame(ped, num.of.err=rep(-1, nr), errate=rep(-1, nr)) 
        tmp <- c()
        for (i in seq_len(nr)) {
            res <- MendErrCal(off.gen.set[i, ], par.gen.set[i, ])  # 调用函数MendErrCal()函数进行检�?
            mend.err.res[i, 3] <- res[[1]]
            mend.err.res[i, 4] <- res[[2]]
            if(res[[2]] <= thre){  # 小于阈值，说明亲子关系记录正确
                tmp <- c(tmp, i)
            }
        } 
        if(length(tmp) == 0){  # 如果tmp长度�?0�? 则提示没有正确的系谱记录
            warning("No right pedigree records!")
        }else{
            ped.right.record <- ped[tmp, ]  # 正确的系谱记�?
            write.table(ped.right.record, "ped_right_record.txt", row.names = FALSE, quote=FALSE)
        }  
        if(length(tmp) == nr){ 
            setwd("../")
            warning("No wrong pedigree records!")
            return()
        }else{
            ped.wrong.record <- mend.err.res[-tmp, ]
            write.table(ped.wrong.record, "ped_wrong_record.txt", row.names = FALSE, quote=FALSE)
            setwd("../")
            if(is.picture){  
                nc <- ncol(gen) - 1
                ErratePlot(mend.err.res[, 4], nc)
                png(file="err_bargraph_01.png", bg="transparent")
                ErratePlot(mend.err.res[, 4], nc)
                dev.off() 
                dev.new()
                ErratePlot(ped.wrong.record[, 4], nc)
                png(file="err_bargraph_02.png", bg="transparent")
                ErratePlot(ped.wrong.record[, 4], nc)
                dev.off()
            }
            return(ped.wrong.record[, 1:2])
        }
    }

    # The following function is designed for the parentage(single) inference
    SinPedInfer <- function(ped, gen, ped.wrong.record, ...){
        ped.right.infer <- ParentInfer(ped, gen, ped.wrong.record)  # 调用ParentInfer()函数，返回推断出的正确系谱记�?
        cat("\n\nRes_02:Pedigree Infer\n")
        if(nrow(ped.right.infer) == 0){  # 如果没有推断出的正确系谱则提�?
            warning("The number of single parent&child pairs inferred correctly is 0!")
            return()
        }else{
            dir.create("ped_infer_result")
            setwd("./ped_infer_result")
            rownames(ped.right.infer) <- NULL
            right.ID.merge.show <- merge(ped.wrong.record, ped.right.infer, by=c("offID"), all.x=TRUE, suffixes = c("W", "R"))  # 通过merge()函数�? 系谱中错误ID和正确ID同时merge展示
            cat("The right pedigree", ngettext(nrow(ped.right.infer), "record", "records"), "by inferred:\n")
            if(is.merge.output){  # if TRUE,则系谱中wrongID && rightID同时write到文件中 
                print(right.ID.merge.show)
                write.table(right.ID.merge.show, "ped_right_infer(merge).txt", quote=FALSE, row.names=FALSE)
            }else{
                print(ped.right.infer)
                write.table(ped.right.infer, "ped_right_infer.txt", quote=FALSE, col.names=FALSE)
            }          
            if(nrow(ped.wrong.record) - nrow(ped.right.infer)){  # if TRUE，则部分错误系谱记录没有推断出正确的亲子关系 
                ped.no.infer <- ped.wrong.record[!ped.wrong.record[, 1] %in% ped.right.infer[, 1], ]
                rownames(ped.no.infer) <- NULL
                write.table(ped.no.infer, "ped_no_infer.txt", quote=FALSE, row.names=FALSE)
                cat("The wrong pedigree", ngettext(d, "record", "records"), "without right inference:\n")
                print(ped.no.infer)
            }
            setwd("../")
        }
    }
                                                    
    # The following function is designed for paternity(biparental) test 
    BipPedCheck <- function(ped, gen, ...){
        GetMendErr <- function(ped.dup, gen, ...){
            mend.err.res <- data.frame(ped.dup, num.of.err=rep(-1, nrow(ped.dup)), errate=rep(-1, nrow(ped.dup)))
            gen.set <- TypeCast(gen[match(as.vector(as.matrix(ped.dup)), gen[, 1]), -1])  # 将系谱中所有个体的基因型按照一定次序集合在一�?
            off.gen <- gen.set[seq_len(i <- nrow(ped.dup)), ]  # 子代的基因型集合�?1：i�?
            sir.gen <- gen.set[(i+1):(2*i), ]  # 父本的基因型集合�?(i+1):(2*i)�?
            dam.gen <- gen.set[(2*i+1):(3*i), ]  # 母本的基因型集合�?(2*i+1):(3*i)�?
            if(i==1){
                res <- MendErrCal(off.gen, sir.gen, dam.gen)
                mend.err.res[, 4] <- res$num.of.err
                mend.err.res[, 5] <- res$errate                
            }else{
                for(j in seq_len(i)){
                    res <- MendErrCal(off.gen[j, ], sir.gen[j, ], dam.gen[j, ])
                    mend.err.res[j, 4] <- res$num.of.err
                    mend.err.res[j, 5] <- res$errate
                }
            }
            mend.err.res[, 5] <- sprintf("%.3f", mend.err.res[, 5])
            mend.err.res <- mend.err.res[order(as.numeric(mend.err.res[, 1])), ]
            return(mend.err.res)
        }
        cat("Res_02:Pedigree Check\n")
        if(any(offID.rep.index <- duplicated(ped[, 1]))){  # 判断系谱中是否有offID重复的情�?                  
            offID.rep.record <- ped[ped[, 1] %in% ped[offID.rep.index, 1], ]  # offID重复的记�?
            cat("Warnings:", ngettext(length(rep.offID <- ped[offID.rep.index, 1]), "There is an offID duplicated", "There are some offIDs duplicated"), "in the pedigree:", paste(sQuote(rep.offID), collapse = ","), "\n")
            if(any(tmp <- duplicated(offID.rep.record))){  # 检查系谱中是否有重复记�?
                cat("Moreover, there are even some duplicated records among them:\n")
                print(offID.rep.record[offID.rep.record[, 1] %in% offID.rep.record[tmp, 1], ])
            }
            mend.err.res <- GetMendErr(offID.rep.record, gen)
            if(any(tmp <- mend.err.res[, 5] <= thre)){  # if TRUE, 说明包含重复offID的系谱记录中存在正确的系谱记�?
                cat("[1] The right pedigree ", ngettext(sum(tmp), "record", "records")," with the duplicated offID as follows:\n")
                print(mend.err.res[tmp, ])
                if(sum(!tmp)){
                    cat("[2] The wrong pedigree ", ngettext(sum(!tmp), "record", "records"), " with the duplicated offID as follows:\n")
                    print(mend.err.res[!tmp, ])
                }
            }else{
                cat("The wrong pedigree records with duplicated offID as follows:\n")
                print(mend.err.res)
                warning("All the records with duplicated offID are wrong!")
            }
            setwd("../")
            stop("Please check out the pedigree according to the above warning messages!", call.=FALSE)
        }
        if(any(tmp1 <- ped[, 2] %in% ped[, 3])){  # 判断是否有ID号同时出现在系谱sirID和damID中的情况，if so 必定存在错误系谱记录
            out.of.ped <- ped[tmp1 | ped[, 3] %in% ped[, 2], ]  # 包含上述情况ID号的系谱记录
            mixed.id <- unique(ped[tmp1, 2])
            message("Warning: ", ngettext(length(mixed.id), "There is an ID",  "There are some IDs"), " mixed in both sirID and damID\n",   
                "[1] The total number of rows involving in the mixed ID: ", sQuote(nrow(out.of.ped)), "\n\n", 
                "[2] The mixed ID recommmended to check: ", paste(sQuote(mixed.id), collapse=","), "\n\n",
                "[3] The rows recommmended to check:"
                )
            if(nrow(out.of.ped[out.of.ped[, 2] %in% mixed.id, ]) <= nrow(out.of.ped[out.of.ped[, 3] %in% mixed.id, ])){  # 推荐需要进行校正的系谱记录�?
                print(out.of.ped[out.of.ped[, 2] %in% mixed.id, ]) 
            }else{
                print(out.of.ped[out.of.ped[, 3] %in% mixed.id, ])  
            }  
            setwd("../")
            stop("Please check out the pedigree records according to the above warning messages!", call.=FALSE)
        } 
        cat("No prompt messages printed\n\n")
        mend.err.res <- GetMendErr(ped, gen)
        write.table(mend.err.res, "mend_error_result.txt", row.names=F, col.names=T, quote=F) 
        ped.wrong.record <- mend.err.res[mend.err.res[, 5] > thre, ]  
        dir.create("ped_check_result", showWarnings=TRUE)   
        setwd("./ped_check_result")  
        ped.right.record <- mend.err.res[mend.err.res[, 5] <= thre, ]       
        if(nrow(ped.right.record) == 0){  
            warning("No right pedigree records!")
        }else{
            write.table(ped.right.record, "ped_right_record.txt", row.names = FALSE, quote=FALSE)
        }  
        if(nrow(ped.wrong.record) == 0){ 
            setwd("../")
            warning("No wrong pedigree records!")
            return()   
        }else{
            ped.wrong.record <- ped.wrong.record[order(as.numeric(ped.wrong.record[, 1])), ] 
            write.table(ped.wrong.record, "ped_wrong_record.txt", row.names = FALSE, quote=FALSE)
            setwd("../")
            if(is.picture){
                nc <- ncol(gen)-1    
                ErratePlot(mend.err.res[, 5], nc)
                png(file="err_bargraph_01.png", bg="transparent")
                ErratePlot(mend.err.res[, 5], nc)
                dev.off()     
                dev.new()
                ErratePlot(ped.wrong.record[, 5], nc)  
                png(file="err_bargraph_02.png", bg="transparent")
                ErratePlot(ped.wrong.record[, 5], nc) 
                dev.off()
            }
            return(ped.wrong.record[, 1:3])
        } 
    }

    # The following function is designed for the parentage(biparental) inference
    BipPedInfer <- function(ped, gen, ped.wrong.record, ...){
        sir.infer.res <- ParentInfer(ped, gen, ped.wrong.record)  # 调用ParentInfer()函数, 返回值为推断出的正确系谱记录(offspring && sire)
        dam.infer.res <- ParentInfer(ped, gen, ped.wrong.record, ped.col=3)  # 同上分析(返回offspring && dam)
        ped.right.infer <- merge(sir.infer.res$tmp2, dam.infer.res$tmp2, by=c("offID"), sort=TRUE)  
        damID.wrong.record <- sir.infer.res$tmp3  
        sirID.wrong.record <- dam.infer.res$tmp3 
        dir.create("ped_infer_result", showWarnings=TRUE) 
        setwd("./ped_infer_result")
        cat("Res_03:Pedigree Infer\n")
        if(nrow(ped.right.infer) == 0){  # 如果没有推断正确的系谱记�?
            if(nd <- nrow(damID.wrong.record)){  
                write.table(damID.wrong.record, "damID_wrong_record.txt", row.names=FALSE, quote=FALSE)
                cat(ngettext(nd, "\nA record with damID wrong:\n", "\nRecords with damID wrong:\n"))
                print(damID.wrong.record)
            }
            if(ns <- nrow(sirID.wrong.record)){  
                write.table(sirID.wrong.record, "sirID_wrong_record.txt", row.names=FALSE, qute=FALSE)
                cat(ngettext(ns, "\nA record with sirID wrong:\n", "\nRecords with sirID wrong:\n"))
                print(sirID.wrong.record)
            }
            ped.right.record <- ped[!ped[, 1] %in% ped.wrong.record[, 1], ]
            offID.wrong.record <- merge(ped.wrong.record, ped.right.record, by=c("sirID", "damID"))
            offID.wrong.record <- offID.wrong.record[, c(3, 1, 2)]
            colnames(offID.wrong.record) <- NULL
            if(no <- nrow(offID.wrong.record)){
                write.table(offID.wrong.record, "offID_wrong_record.txt", row.names=FALSE, qute=FALSE)
                cat(ngettext(ns, "\nA record with offID wrong:\n", "\nRecords with offID wrong:\n"))
                print(offID.wrong.record)
            }
            if(np <- nrow(ped.wrong.record) - sum(nd, ns, no)){ 
                parIDs.all.wrong.record <- ped.wrong.record[!ped.wrong.record[, 1] %in% c(damID.wrong.record[, 1], sirID.wrong.record[, 1], offID.wrong.record[, 1]), ]
                write.table(allID.wrong.record, "parIDs_all_wrong_record.txt", row.names=FALSE, quote=FALSE) 
                cat(ngettext(np, "\nA record with parIDs all wrong:\n", "\nRecords with parIDs all wrong:\n"))
                print(parIDs.all.wrong.record)
            }
            warning("The number of parents && child pairs inferred correctly is 0!")  
        }else{
            nrow.of.ped.wong.record <- nrow(ped.wrong.record)
            if(any(tmp <- ped.right.infer[, 2] %in% ped.right.infer[, 3])){  # 判断ped.right.infer结果中是否存在sirID && damID相同的情�?
                warning(patste0(ngettext(sum(tmp), "\nThere is a  record ", "\nThere are some records"), "in ped.right.infer with the sirID&&damID same!\n"))
                print(ped.right.infer[tmp, ])
                ped.right.infer <- ped.right.infer[!tmp, ]
            }
            write.table(ped.right.infer, "ped_right_infer.txt", row.names = FALSE, quote=FALSE)
            if(any(tmp <- !ped.wrong.record[, 1] %in% ped.right.infer[, 1])){  # 没有找到正确亲子关系的错误系谱记�?
                ped.no.infer <- ped.wrong.record[tmp, ]  
                cat("\nThe wrong", ngettext(sum(tmp), "record", "records"), "without right inference\n")
                print(ped.no.infer)
                write.table(ped.no.infer, "ped_no_infer.txt", row.names = FALSE, quote=FALSE) 
            }
            if(nrow(damID.wrong.record)){
                right.damID.merge.show <- merge(damID.wrong.record, ped.right.infer, by=c("offID", "sirID"), all.x=TRUE, suffixes = c("W", "R"))  # 通过merge()函数将正确的 damID和错误的damID合并展示
                if(any(is.na(right.damID.merge.show[, 4]))){  # 如果存在NA的情况，则将包含NA的行排在后面
                    right.damID.merge.show <- right.damID.merge.show[order(as.numeric(right.damID.merge.show[, 4])), ]              
                }
                cat(paste("\nThe only_damID_wrong", ngettext(nrow(damID.wrong.record), "record:", "records:")), "\n")
                print(right.damID.merge.show)
                if(is.merge.output){
                    write.table(right.damID.merge.show, "damID_wrong_record.txt", quote=FALSE, row.names=FALSE) 
                }else{
                    write.table(damID.wrong.record, "damID_wrong_record.txt", quote=FALSE, row.names=FALSE)                         
                }
                ped.right.infer <- ped.right.infer[!ped.right.infer[, 1] %in% damID.wrong.record[, 1], ]  # 为避免ped.right.infer中同一记录多次参与比较，所以将上述参与damID.wrong.record比较的记录删�?                
                ped.wrong.record <- ped.wrong.record[!ped.wrong.record[, 1] %in% damID.wrong.record[, 1], ]
            }
            if(nrow(sirID.wrong.record)){
                right.sirID.merge.show <- merge(sirID.wrong.record, ped.right.infer, by=c("offID", "damID"), all.x=TRUE, suffixes = c("W", "R"))  # 同上分析
                if(any(is.na(right.sirID.merge.show[, 4]))){  # 同上分析
                    right.sirID.merge.show <- right.sirID.merge.show[order(as.numeric(right.sirID.merge.show[, 4])), ]              
                }          
                cat(paste("\nThe only_sirID_wrong", ngettext(nrow(sirID.wrong.record), "record:", "records:")), "\n")
                print(right.sirID.merge.show)
                if(is.merge.output){
                    write.table(right.sirID.merge.show, "sirID_wrong_record.txt", quote=FALSE, row.names=FALSE) 
                }else{
                    write.table(sirID.wrong.record, "sirID_wrong_record.txt", quote=FALSE, row.names=FALSE)                         
                }
                ped.right.infer <- ped.right.infer[!ped.right.infer[, 1] %in% sirID.wrong.record[, 1], ]  
                ped.wrong.record <- ped.wrong.record[!ped.wrong.record[, 1] %in% sirID.wrong.record[, 1], ]
            }
            if(nrow(ped.wrong.record)){
                right.offID.merge.show <- merge(ped.wrong.record, ped.right.infer, by=c("sirID", "damID"), suffixes = c("W", "R"))  
                offID.wrong.record <- right.offID.merge.show[, c(3, 1, 2)]
                if(nrow(right.offID.merge.show)){
                    colnames(offID.wrong.record) <- c("offID", "sirID", "damID")
                    cat(paste("\nThe only_offID_wrong", ngettext(nrow(offID.wrong.record), "record:", "records:")), "\n")
                    print(right.offID.merge.show)
                    if(is.merge.output){
                        write.table(right.offID.merge.show, "offID_wrong_record.txt", quote=FALSE, row.names=FALSE) 
                    }else{
                        write.table(offID.wrong.record, "offID_wrong_record.txt", quote=FALSE, row.names=FALSE)                         
                    }     
                    ped.right.infer <- ped.right.infer[!ped.right.infer[, 1] %in% offID.wrong.record[, 1], ]
                    ped.wrong.record <- ped.wrong.record[!ped.wrong.record[, 1] %in% offID.wrong.record[, 1], ]                    
                }
            }
            if(nrow(ped.wrong.record)){
                parIDs.all.wrong.record <- ped.wrong.record[!ped.wrong.record[, 1] %in% offID.wrong.record[, 1], ]  # 系谱亲本ID号都记录错误的系谱记�?
                rownames(parIDs.all.wrong.record) <- NULL
                right.parIDs.merge.show <- merge(parIDs.all.wrong.record, ped.right.infer, by=c("offID"), suffixes = c("W", "R"))
                cat("\nThe parIDs_all_wrong", ngettext(nrow(parIDs.all.wrong.record), "record:\n", "records:\n"))
                print(right.parIDs.merge.show)
                if(is.merge.output){
                    write.table(right.parIDs.merge.show, "parIDs_all_wrong_record.txt", quote=FALSE, row.names=FALSE) 
                }else{
                    write.table(parIDs.all.wrong.record, "parIDs_all_wrong_record.txt", quote=FALSE, row.names=FALSE)                         
                }
            }
            cat("The total: ", nrow.of.ped.wong.record, "\n\n")  
        }
        setwd("../")
    }

    ##============================ The program executes here ===============================
    options(digits = 3)
    ped.name <- colnames(ped)
    if(is.null(ped.name)){  # 检验输入的系谱有没有列�?
        stop("The colnames of the pedigree are missing!")
    }else if(paste(sQuote(ped.name), collapse=",") != (name.tmp <- ngettext(ncol(ped)-1, paste(sQuote(c("offID", "parID")), collapse=","), paste(sQuote(c("offID","sirID","damID")), collapse=",")))) {  # 系谱按照能且只能按照"offID"�?"sirID", "damID"三列的形式依次排�?
        stop("The colnames of the pedigree should be corresponding to the order of following name:", name.tmp)
    }
    gen.symbol <- unique(as.vector(as.matrix(gen[, -1])))
    if(all(gen.symbol %in% c(1, 2, 3, 5))){  # 表示SNP位点基因的符号为'1', '2', '3', '5'分别表示：aa, Aa, AA, unknown; 否则出错
        mode <- ncol(ped)-1
        ped <- SnpDataCheck(ped, gen)
        if(mode == 1){  # 'mode==1'表示单亲情况下的亲子鉴定
            dir.create("paternity_test(Single)", showWarnings=TRUE) 
            setwd("./paternity_test(Single)")
            ped.wrong.record <- SinPedCheck(ped, gen) 
            if(!is.null(ped.wrong.record)){
                SinPedInfer(ped, gen, ped.wrong.record)  
            }  
            setwd("../")  
        }else if(mode == 2){  # 'mode==2'表示双亲情况下的亲子鉴定
            dir.create("paternity_test(Biparental)", showWarnings=TRUE)
            setwd("./paternity_test(Biparental)")
            ped.wrong.record <- BipPedCheck(ped, gen)   
            if(!is.null(ped.wrong.record)){
                BipPedInfer(ped, gen, ped.wrong.record)  
            }
            setwd("../")  
        }else{
            stop("Please select the right mode: '1' representing 'single' or '2' representing 'biparental'!")
        }             
    }else{
        stop("The symbol of genotype is wrong: only '1', '2', '3' and '5' avaliable!")
    }
}  

##============================ program starts here ==================
setwd("C:/")            # set your work directory
gen <- as.matrix(read.table("gen.txt", header=F))   # read genotype file
ped <- as.matrix(read.table("ped.txt", header=T))   # read pedigree file
ParTest(ped, gen, thre...)  # paternity test, thre= threshold value.

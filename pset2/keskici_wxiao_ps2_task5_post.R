#! /usr/bin/env Rscript

source('keskici_wxiao_ps2_functions.R')

TASK.NUM = 5
theta.draws <- 3    # theta.nsims * N is the # of total simulations we'll run
pair.nums <- 4

aggregate.cover(pair.nums, theta.draws, TASK.NUM, TRUE)

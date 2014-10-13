#! /usr/bin/env Rscript

source('keskici_wxiao_ps2_functions.R')

TASK.NUM = 3

theta.draws <- 9
pair.nums <- 4

aggregate.cover(pair.nums, theta.draws, TASK.NUM)

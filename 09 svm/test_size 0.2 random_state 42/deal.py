#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :deal.py
# @Time :2022/2/10 15:04

deal_file = open('train suoyin.txt', "r")
zong_file = open('all sample.txt', 'r')
new_file = open('train sample.txt', 'w')

zong_dish = {}
for line in zong_file.readlines():
    id = line.split()[0].split()[0]
    num = line.split()[1:]
    zong_dish[str(id)] = num

for line in deal_file.readlines():
    new_line = line.split('\n')[0]
    if str(new_line) in zong_dish.keys() and zong_dish[new_line] != []:
        new_file.write(new_line + "\t" + "\t".join(zong_dish[new_line]) + '\n')


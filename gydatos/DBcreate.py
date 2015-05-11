#!/usr/bin/env Python

import sqlite3 as sq

cnx = sq.connect('hydro_whales.db')

print ('Connector:', cnx)

cnx.execute("DROP TABLE NEXTID")
cnx.execute("DROP TABLE TRIAD")
cnx.execute("DROP TABLE WHALEDETECTION")

cnx.execute(" CREATE TABLE NEXTID ( \
IDName varchar(255),\
current int)")

cnx.execute(" CREATE TABLE WHALEDETECTION ( \
detid int, \
StaName varchar(255),\
time real,\
time_cc  real,\
whale_disc int)")


cnx.execute(" CREATE TABLE TRIAD ( \
TriadName varchar(255),\
triadid int,\
detid1 int, \
detid2 int, \
detid3 int, \
AbsDirection real,\
Absfit  real,\
CCDirection real,\
CCfit  real,\
StkDirection real,\
Stkfit  real\
)")

cnx.close()

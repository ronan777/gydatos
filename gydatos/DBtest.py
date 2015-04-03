#!/usr/bin/env Python

import _mysql

cnx = _mysql.connect(host='sql4.freesqldatabase.com',user='sql457121', passwd='fN2%aL4!', db='sql457121')

print ('Connector:', cnx)

cnx.query("""
CREATE TABLE Persons
(
PersonID int,
LastName varchar(255),
FirstName varchar(255),
Address varchar(255),
City varchar(255))""")

cnx.close()

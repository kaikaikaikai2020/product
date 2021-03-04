# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 22:11:13 2021

@author: Asus
"""

import pandas as pd
import glob
import os

list_of_files = glob.glob('L:\Dropbox\Dropbox\project folder from my asua computer\Project\S60program rl pair\program\S60P3_para\position*-800.xlsx')
'''
print (list_of_files)
'''

import datetime
from pandas.tseries.offsets import BDay

today = datetime.datetime.today()
today_str = today.strftime('%Y-%m-%d')
yst = (today-BDay(1))
yst_str = yst.strftime('%Y-%m-%d')
final = pd.DataFrame()
for file in list_of_files:
    df = pd.read_excel(file)
    df['filename']= file
    #df = df.loc[df['trade_time']==yst_str|df['trade_time']==today_str]
    df1 = df.loc[df['trade_time']==yst_str]
    final = final.append(df1)
    df2 = df.loc[df['trade_time']==today_str]
    final = final.append(df2)
print(final)
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from smtplib import SMTP
import smtplib
import sys

reci = ['caifengsteven@gmail.com']
emaillist = [elem.strip().split(',') for elem in reci] 

msg = MIMEMultipart()
msg['Subject'] = "daily pair trading " +today_str
msg['From']='caifengsteven@gmail.com'

html = """\
<html>
  <head></head>
  <body>
    {0}
  </body>
</html>
""".format(final.to_html())

part1 = MIMEText(html, 'html')
msg.attach(part1)

server = smtplib.SMTP('smtp.gmail.com', 587)
server.starttls()
server.login(msg['From'], "352471Cf")
server.sendmail(msg['From'], emaillist , msg.as_string())





   
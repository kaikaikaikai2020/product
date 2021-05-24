# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:32:57 2021

@author: ASUS
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 08:31:50 2020

@author: ASUS
"""
from git import Repo
import pandas as pd
import glob
import os

#Repo.clone_from('https://github.com/kaikaikaikai2020/test.git','g:/adr')
'''from bs4 import BeautifulSoup
import requests
x=requests.get('https://github.com/kaikaikaikai2020/test/tree/master/ADR/ADR/ZX02-program/program/%E8%AE%A1%E7%AE%97%E7%BB%93%E6%9E%9C')
soup = BeautifulSoup(x.content, 'html.parser')
#print (soup.find_all('a'))
result = [ xls['href'] for xls in soup.find_all('a') if 'csv' in xls['href']]
for url in result:
    url = 'http://github.com'+url
    #print(url)
    downloaded_file = requests.get(url)
    print(downloaded_file)
    if os.path.isfile(downloaded_file):
        print ("File exist")
    else:
        print ("File not exist")
#print(result[0])
'''

#list_of_files = glob.glob('G:/adr/ADR/ADR/ZX02-program/program/计算结果/*.csv') # * means all if need specific format then *.csv
list_of_files = glob.glob('G:/github repo/product/ADR/ADR/ZX02-program/program/计算结果/支线任务2-*操作明细-1.50.csv')
#list_of_files = glob.glob('G:\计算结果\支线任务2-*操作明细-1.50.csv')
latest_file = max(list_of_files, key=os.path.getctime)
print(latest_file)
df = pd.read_csv(latest_file, dtype={'date':'string'})
print(df)
df['multiple']=1
#df = df.set_index('date')

import datetime
from pandas.tseries.offsets import BDay

today = datetime.datetime.today()
#need to change the date 

today_str= today.strftime('%Y-%m-%d')


'''

calculate the level to buy and sell those ADR/DR names

'''

df = df.loc[df['date'] == today_str]
df = df.loc[df['methodID'] == 'signal1-model2']
df.loc[df['pairName'] == 'INFO-INFY',['multiple']] =1
df.loc[df['pairName'] == '2330-TSM',['multiple']] =5
df.loc[df['pairName'] == 'ICICIBC-IBN',['multiple']] =2
df.loc[df['pairName'] == 'TTMT-TTM',['multiple']] =5
df.loc[df['pairName'] == 'HDFCB-HDB',['multiple']] =3
df.loc[df['pairName'] == '2303-UMC',['multiple']] =5

df['Sell_US_level'] = df['asiaOP']/df['exOP']*df['multiple']/(1-df['HighBond'])
df['Buy_US_level'] = df['asiaOP']/df['exOP']*df['multiple']/(1-df['LowBond'])
print(df)
df1 = df.loc[(df['pairName']=='2303-UMC')|(df['pairName']=='2330-TSM')]
df1= df1[['pairName','date', 'Buy_US_level','Sell_US_level']]

df2 = df.loc[(df['pairName']=='INFO-INFY')|(df['pairName']=='ICICIBC-IBN')|(df['pairName']=='TTMT-TTM')|(df['pairName']=='HDFCB-HDB')]
df2= df2[['pairName','date', 'Buy_US_level','Sell_US_level']]


from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from smtplib import SMTP
import smtplib
import sys


recipients = ['caifengsteven@gmail.com','steven.cai@nomura.com'] 
emaillist = [elem.strip().split(',') for elem in recipients]
msg = MIMEMultipart()
msg['Subject'] = "daily ADR levels US open Taiwan"+today_str
msg['From'] = 'caifengsteven@gmail.com'


html = """\
<html>
  <head></head>
  <body>
    {0}
  </body>
</html>
""".format(df1.to_html())

part1 = MIMEText(html, 'html')
msg.attach(part1)

server = smtplib.SMTP('smtp.gmail.com', 587)
server.starttls()
server.login(msg['From'], "352471Cf")
server.sendmail(msg['From'], emaillist , msg.as_string())


recipients = ['caifengsteven@gmail.com','steven.cai@nomura.com'] 
emaillist = [elem.strip().split(',') for elem in recipients]
msg = MIMEMultipart()
msg['Subject'] = "daily ADR levels US open India"+today_str
msg['From'] = 'caifengsteven@gmail.com'


html = """\
<html>
  <head></head>
  <body>
    {0}
  </body>
</html>
""".format(df2.to_html())

part1 = MIMEText(html, 'html')
msg.attach(part1)

server = smtplib.SMTP('smtp.gmail.com', 587)
server.starttls()
server.login(msg['From'], "352471Cf")
server.sendmail(msg['From'], emaillist , msg.as_string())


'''

sending out email

'''


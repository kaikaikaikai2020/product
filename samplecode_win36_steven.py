# -*- coding: utf-8 -*-
from dataapi_win36 import Client
if __name__ == "__main__":
    try:
        client = Client()
        client.init('473c6c98362ce267eafb83a9bf37bcd64eb25d739829b1356c426bc1d3f35ac7')
        from datetime import timedelta, date

        def daterange(start_date, end_date):
            for n in range(int((end_date - start_date).days)):
                yield start_date + timedelta(n)

        start_date = date(2018, 1, 1)
        end_date = date(2020, 12, 10)
        for single_date in daterange(start_date, end_date):
            print(single_date.strftime("%Y%m%d"))
            current_date = single_date.strftime("%Y%m%d")
		#调用股票基本信息(仅供参考)
            #s = '{}_{:05}_{}_{:.5f}'.format(s1, i, s2, f)
            #url1='/api/market/getHKMktTuoHistOneDay.json?dataDate={}&tradeDirection=NB'.format(current_date) 
            #print(url1)
            #url1='/api/market/getHKMktTuoHistOneDay.json?dataDate=20190524&tradeDirection=NB'
            #code, result = client.getData(url1)#调用getData函数获取数据，数据以字符串的形式返回
            #if code==200:
            #    print(result.decode('utf-8'))#url1须为json格式，才可使用utf-8编码
            #    filename = "getSec{}.csv".format(current_date)
            #    file_object = open(filename, 'w')
            #    file_object.write(result.decode('GB18030'))#url2须为csv格式才可使用GB18030编码，写入到getSecST.csv文件
            #    file_object.close( )
    			#pd_data=pd.DataFrame(eval(result)['data'])#将数据转化为DataFrame格式
            #else:
            #    print (code)
            #    print (result)
            
            #url2='/api/market/getHKMktTuoHistOneDay.csv?dataDate=20190524&tradeDirection=NB'
            url2='/api/market/getHKMktTuoHistOneDay.csv?dataDate={}&tradeDirection=NB'.format(current_date) 
            code, result = client.getData(url2)
            if(code==200):
                filename = "getSec{}.csv".format(current_date)
                file_object = open(filename, 'w')
                #file_object = open('getSecST.csv', 'w')
                file_object.write(result.decode('GB18030'))#url2须为csv格式才可使用GB18030编码，写入到getSecST.csv文件
                file_object.close( )
            else:
                print (code)
                print (result) 
    except Exception as e:
        #traceback.print_exc()
        raise e
    
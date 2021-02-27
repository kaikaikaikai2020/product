# -*- coding: utf-8 -*-
from dataapi_win36 import Client
if __name__ == "__main__":
    try:
        client = Client()
        client.init('473c6c98362ce267eafb83a9bf37bcd64eb25d739829b1356c426bc1d3f35ac7')
		#调用股票基本信息(仅供参考)
        url1='/api/market/getHKMktTuoHistOneDay.json?dataDate=20190524&tradeDirection=NB'
        code, result = client.getData(url1)#调用getData函数获取数据，数据以字符串的形式返回
        if code==200:
            print(result.decode('utf-8'))#url1须为json格式，才可使用utf-8编码
			#pd_data=pd.DataFrame(eval(result)['data'])#将数据转化为DataFrame格式
        else:
            print (code)
            print (result)
		#调取沪深股票ST标记数据
        url2='/api/market/getHKMktTuoHistOneDay.csv?dataDate=20190524&tradeDirection=NB'
        code, result = client.getData(url2)
        if(code==200):
            file_object = open('getSecST.csv', 'w')
            file_object.write(result.decode('GB18030'))#url2须为csv格式才可使用GB18030编码，写入到getSecST.csv文件
            file_object.close( )
        else:
            print (code)
            print (result) 
    except Exception as e:
        #traceback.print_exc()
        raise e
    
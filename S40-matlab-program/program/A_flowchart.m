clear
%将预测者分钟数据转换为时间序列数据（后续可以考虑表分区，一张表就可以）
title_str = '将预测者分钟数据转换为时间序列数据（后续可以考虑表分区，一张表就可以）';
order_str = 'M_YCZmin_section_to_series';
run_program_adair(order_str,title_str);

%A股指数 已经对接
title_str =  'A股指数 已经对接';
order_str = 'M_Astock_Indicator';
run_program_adair(order_str,title_str);

%A股商品期货
title_str =  'A股商品期货';
order_str = 'M_future_CF';
run_program_adair(order_str,title_str);

%股指期货
title_str =  '股指期货';
order_str = 'M_Astock_Future_Indicator';
run_program_adair(order_str,title_str);

%外汇
title_str =  '外汇';
order_str = 'M_Foreign_exchange';
run_program_adair(order_str,title_str);
%A股及指数 我的硬盘快不行了，正常的硬盘16核心需要13分钟，8核心应该再26-30分钟能算完。
%下了功夫优化速度，原来需要60分钟。
title_str =  'A股及指数';
order_str ='M_Astock_final';
run_program_adair(order_str,title_str);

%世界主要指数数据 自动判断交易时间
title_str =  '世界主要指数数据 自动判断交易时间';
order_str ='M_Index_Foreign_main_index_update2';
run_program_adair(order_str,title_str);

%指数测试加权
title_str =  '指数测试加权';
order_str ='M_Astock_Indicator_weightupdate';
run_program_adair(order_str,title_str);

%美股
title_str =  '美股';
order_str ='M_Astock_update2_american';
run_program_adair(order_str,title_str);

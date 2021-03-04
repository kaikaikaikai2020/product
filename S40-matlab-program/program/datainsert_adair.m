function OK = datainsert_adair(tablename,colnames,data)
    conn = mysql_conn();
    try
        datainsert(conn,tablename,colnames,data);
        OK = true;
    catch
        sprintf('datainsert error!')
        OK = false;
    end
    close(conn)
end
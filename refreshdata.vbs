Option Explicit



Dim xlApp, xlBook, xlSheet



Set xlApp = CreateObject("Excel.Application")

Set xlBook = xlApp.Workbooks.Open("C:\Users\ASUS\Documents\GitHub\test\Non Fungible ADR strategies.xlsx")

Set xlSheet = xlBook.worksheets.item(1)

xlBook.RefreshAll
xlBook.Save

xlBook.Close

xlApp.Quit



Set xlSheet = Nothing

Set xlBook = Nothing

Set xlApp = Nothing
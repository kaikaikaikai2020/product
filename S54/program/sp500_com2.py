# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 22:40:34 2020

@author: adair2019
"""

__author__ = "Curtis Miller"
__copyright__ = "Copyright (c) 2017, Curtis Grant Miller"
__credits__ = ["Curtis Miller"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Curtis Miller"
__email__ = "cgmil@msn.com"
__status__ = "Experimental"
 
import pandas as pd
from pandas import DataFrame
import argparse
import quandl
import pandas_datareader as web
from time import sleep
import datetime as dt
import sys
 
def get_sp500_data(start=dt.datetime.strptime("1997-01-01", "%Y-%m-%d"),
                   end=dt.datetime.now(), use_quandl=True, adjust=True, inner=True,
                   sleeptime=2, verbose=True):
    """Fetches S&P 500 data
     
    args:
        start: datetime; The earliest possible date
        end: datetime; The last possible date
        use_quandl: bool; Whether to fetch data from Quandl (reverts to Google if False)
        adjust: bool; Whether to use adjusted close (only works with Quandl)
        inner: bool; Whether to use an inner join or outer join when combining series (inner has no missing data)
        sleeptime: int; How long to sleep between fetches
        verbose: bool; Whether to print a log while fetching data
     
    return:
        DataFrame: Contains stock price data
    """
     
    join = "outer"
    if inner:
        join = "inner"
     
    symbols_table = pd.read_html("https://en.wikipedia.org/wiki/List_of_S%26P_500_companies",
                                 header=0)[0]
    symbols = list(symbols_table.loc[:, "Ticker symbol"]) # 可能是页面变化，这里改成Symbol
 
    sp500 = None
    for s in symbols:
        sleep(sleeptime)
        if verbose:
            print("Processing: " + s + "...", end='')
        try:
            if use_quandl:
                s_data = quandl.get("WIKI/" + s, start_date=start, end_date=end)
                if adjust:
                    s_data = s_data.loc[:, "Adj. Close"]
                else:
                    s_data = s_data.loc[:, "Close"]
            else:
                s_data = web.DataReader(s, "google", start, end).loc[:, "Close"]
            s_data.name = s
            s_data.dropna()
            if s_data.shape[0] > 1:
                if sp500 is None:
                    sp500 = DataFrame(s_data)
                else:
                    sp500 = sp500.join(s_data, how=join)
                if verbose:
                    print(" Got it! From", s_data.index[0], "to", s_data.index[-1])
            else:
                if verbose:
                    print(" Sorry, but not this one!")
        except Exception:
            if verbose:
                print(" Sorry, but not this one!")
 
    badsymbols = list(set(s) - set(sp500.columns))
    if verbose and len(badsymbols) > 0:
        print("There were", len(badsymbols), "symbols for which data could not be obtained.")
        print("They are:", ", ".join(badsymbols))
     
    return sp500
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fetches S&P 500 data")
    parser.add_argument("-v", "--verbose", action="store_true", default=True, dest="verbose",
                        help="Print extra output [default]")
    parser.add_argument("--quietly", action="store_false",
                        dest="verbose", help="Don't print extra output")
    parser.add_argument("-f", "--file", type=str, dest="csv_name",
                        default="sp-500.csv",
                        help="CSV file to save data to [default: sp-500.csv]")
    parser.add_argument("-s", "--sleep", type=int, dest="sleeptime",
                        default=2,
                        help="Time (seconds) between fetching symbols [default: 2] (don't flood websites with requests!)")
    parser.add_argument("--inner", action="store_true", default=False, dest="inner",
                        help="Inner join; only dates where all symbols have data will be included")
    parser.add_argument("--start", type=str, dest="start",
                        default="1997-01-01",
                        help="Earliest date (YYYY-MM-DD) to include [default: 1997-01-01]")
    parser.add_argument("--end", type=str, dest="end",
                        default="today",
                        help='Last date (YYYY-MM-DD or "today") to include [default: "today"]')
    # parser.add_argument("-k", "--key", type="character", dest="api_key",
    #                     default=NULL,
    #                     help="Quandl API key, needed if getting Quandl data")
    parser.add_argument("-q", "--quandl", action="store_true", default=False,
                        dest="use_quandl", help="Get data from Quandl")
    parser.add_argument("-a", "--adjust", action="store_true", default=False,
                        dest="adjust", help="Adjust prices (Quandl only)")
    parser.add_argument("--about", action="store_true", default=False,
                        dest="about",
                        help="Print information about the script and its usage, then quit")
 
    args = parser.parse_args()
 
    if args.about:
        print(sys.argv[0], "\n(c) 2017 Curtis Miller\n",
          "Licensed under GNU GPL v. 3.0 available at ",
          "https://www.gnu.org/licenses/gpl-3.0.en.html \n",
          "E-mail: cgmil@msn.com\n\n",
          "This script fetches closing price data for ticker symbols included",
          "in the S&P 500 stock index. A list of symbols included in the index",
          "is fetched from this webpage:",
          "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies  The list",
          "is parsed and the symbols included in the list are fetched from",
          "either Google Finance (the default) or Quandl.",
          "If Quandl is the data source, adjusted data can be",
          "fetched instead. The resulting data set is then saved to a CSV",
          "file in the current working directory.\n\n",
          "This package requires the following Python packages be installed in order",
          "to work (all of which are available through pip):\n\n",
          "* pandas\n",
          "* pandas-datareader\n",
          "* quandl\n\n",
          "This script was written by Curtis Miller and was made available on ",
          "his website: https://ntguardian.wordpress.com\n\n",
          "You can read more about this script in the following article: ",
          "https://ntguardian.wordpress.com/blog\n\n")
        quit()
 
    if args.end == "today":
        args.end = dt.datetime.now()
    else:
        args.end = dt.datetime.strptime(args.end, "%Y-%m-%d")
    args.start = dt.datetime.strptime(args.start, "%Y-%m-%d")
    sp500 = get_sp500_data(start=args.start,
                   end=args.end, use_quandl=args.use_quandl, adjust=args.adjust, inner=args.inner,
                   sleeptime=args.sleeptime, verbose=args.verbose)
    sp500.to_csv(args.csv_name)
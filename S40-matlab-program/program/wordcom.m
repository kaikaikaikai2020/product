%class for word com in Matlab
%methods included:
%create word
%insert table
%insert figure
%insert strings
classdef wordcom < handle
    properties
        WH        % WordHandel
        AWH       % ActWordhandle
        PositionS %
        FNs       % word filenames
        Tag = 'wordcom';
        fnPDF
    end
    methods
        %打开word
        function obj = wordcom(word_file_p)
            % Start an ActiveX session with Word:
            actx_word = actxserver('Word.Application');
            %actx_word.Visible = false;
            actx_word.Visible = true;
            trace(actx_word.Visible);
            if exist(word_file_p,'file');
            	delete(word_file_p);
            end
            if ~exist(word_file_p,'file');  % Create new document:
                word_handle = invoke(actx_word.Documents,'Add');
            else  % Open existing document:
                word_handle = invoke(actx_word.Documents,'Open',word_file_p);
            end
            obj.AWH = actx_word;
            obj.WH = word_handle;
            obj.FNs = word_file_p;
        end
        %设定页面
        function  setWordPadge(obj,t,b,l,r)
            if nargin < 5
                t = 2.5;b = 2.5;l = 3.2;r = 3.2;
            end
            value1 = 28.3465;
            obj.WH.PageSetup.TopMargin = t*value1;
            obj.WH.PageSetup.BottomMargin = b*value1;
            obj.WH.PageSetup.LeftMargin = l*value1;
            obj.WH.PageSetup.RightMargin = r*value1;
        end
        %写入文字
        function WordText(obj,inputstr,stylesel)
            actx_word_p = obj.AWH;
            text_p = inputstr;
            if eq(stylesel,1)
                style_p = '标题';
            else
                style_p = '正文';
            end
            enters_p = [0,1];
            for i = 1:enters_p(1)
                actx_word_p.Selection.TypeParagraph; %enter
            end
            actx_word_p.Selection.Style = style_p;
            if(nargin == 5)%check to see if color_p is defined
                actx_word_p.Selection.Font.Color=color_p;     
            end

            actx_word_p.Selection.TypeText(text_p);
            %actx_word_p.Selection.Font.Color='wdColorAutomatic';%set back to default color
            for k=1:enters_p(2)    
                actx_word_p.Selection.TypeParagraph; %enter
            end
        end
        %设置字体间距
        function setFontsize(obj,value)
            obj.AWH.Selection.Style.Font.Size= value;
        end
        %设置字体
        function setFontName(obj,str)
            obj.AWH.Selection.Style.Font.NameAscii=str;
        end
        %设置间距
        function setSpace(obj,value)
            eval(['obj.WH.Content.ParagraphFormat.Space',num2str(value),';']);
        end
        %插入图片
        function pasteFigure(obj,h,tempstr,tempstrPos)
            if nargin < 4
                tempstrPos = 2;
            end
            if eq(tempstrPos,1)
                lighting phong
                set(h,'Renderer','painters')
                hgexport(h, '-clipboard');
                invoke(obj.AWH.Selection,'Paste');
                obj.AWH.Selection.ParagraphFormat.Alignment='wdAlignParagraphCenter';
                delete(h);
            end
            if nargin>2
                actxword = obj.AWH;
                actxword.Selection.TypeParagraph;
                actxword.Selection.TypeText(tempstr);
                actxword.Selection.ParagraphFormat.Alignment='wdAlignParagraphCenter';
                actxword.Selection.TypeParagraph;
                actxword.Selection.ParagraphFormat.Alignment='wdAlignParagraphLeft';
            end
            if eq(tempstrPos,2)
                lighting phong
                set(h,'Renderer','painters')
                hgexport(h, '-clipboard');
                invoke(obj.AWH.Selection,'Paste');
                obj.AWH.Selection.ParagraphFormat.Alignment='wdAlignParagraphCenter';
                delete(h);
            end
        end
        %插入表格
        function insertTable(obj,tablere,tempstr)
            actxword = obj.AWH;
            if nargin > 2                
                actxword.Selection.TypeParagraph;
                actxword.Selection.TypeText(tempstr);
                actxword.Selection.ParagraphFormat.Alignment='wdAlignParagraphCenter';
            end
            WordCreateTable(actxword,size(tablere,1),size(tablere,2),tablere,1); 
        end
        %转换为pdf
        function transPDF(obj,fnPDF)
            invoke(obj.WH,'ExportAsFixedFormat',fnPDF,'wdExportFormatPDF')
        end
        %关闭word
        function CloseWord(obj)
            actx_word_p = obj.AWH;
            word_handle_p = obj.WH;
            word_file_p = obj.FNs;
            if ~exist(word_file_p,'file')
                % Save file as new:
                invoke(word_handle_p,'SaveAs',word_file_p,1);
            else
                % Save existing file:
                invoke(word_handle_p,'Save');
            end
            % Close the word window:
            invoke(word_handle_p,'Close');            
            % Quit MS Word
            invoke(actx_word_p,'Quit');            
            % Close Word and terminate ActiveX:
            delete(actx_word_p);            
        end
    end
end


function WordCreateTable(actx_word_p,nr_rows_p,nr_cols_p,data_cell_p,enter_p,rvs) 
    %Add a table which auto fits cell's size to contents
    if(enter_p(1))
        actx_word_p.Selection.TypeParagraph; %enter
    end
    %create the table
    %Add = handle Add(handle, handle, int32, int32, Variant(Optional))
    %actx_word_p.ActiveDocument.Tables.Add(actx_word_p.Selection.Range,nr_rows_p,nr_cols_p,0,0);
    actx_word_p.ActiveDocument.Tables.Add(actx_word_p.Selection.Range,nr_rows_p,nr_cols_p);
    %Hard-coded optionals                     
    %first 1 same as DefaultTableBehavior:=wdWord9TableBehavior
    %last  1 same as AutoFitBehavior:= wdAutoFitContent
    DTI = actx_word_p.ActiveDocument.Tables.Item(actx_word_p.ActiveDocument.Tables.Count);
    DTI.Borders.OutsideLineStyle = 'wdLineStyleSingle';
    %DTI.Borders.OutsideLineWidth = 'wdLineWidth150pt';
    DTI.Borders.InsideLineStyle = 'wdLineStyleSingle';
    %DTI.Borders.InsideLineWidth = 'wdLineWidth150pt';
    DTI.Rows.Alignment = 'wdAlignRowCenter';    
    if nargin<6
        rvs = [];
        if nr_rows_p>1
            for i = 1:nr_cols_p
                rvs(i) = length(data_cell_p{2,i});
            end
            if rvs(1)<rvs(2)
                rvs(1)=rvs(2);
            end
        else
            rvs = ones(nr_cols_p,1);
        end    
    end
    
    val = (actx_word_p.activeDocument.pageSetup.PageWidth-actx_word_p.activeDocument.pageSetup.LeftMargin - actx_word_p.activeDocument.pageSetup.RightMargin)/sum(rvs).*rvs*0.8;
    %val(val>=49) = 49;
    val(:) = 20;
    for i = 1:nr_cols_p;DTI.Columns.Item(i).Width = val(i);end
    %write the data into the table
    for r=1:nr_rows_p
        for c=1:nr_cols_p
            %write data into current cell
            if eq(c,1)
                WordText0(actx_word_p,data_cell_p{r,c},'Normal',[0,0],'wdAlignParagraphLeft');
                
            else
                WordText0(actx_word_p,data_cell_p{r,c},'Normal',[0,0],'wdAlignParagraphRight');
            end
            if(r*c==nr_rows_p*nr_cols_p)
                %we are done, leave the table
                actx_word_p.Selection.MoveDown;
            else%move on to next cell 
                actx_word_p.Selection.MoveRight;
            end            
        end
    end
end

function WordText0(actx_word_p,text_p,style_p,enters_p,sel0,color_p)
%if ~isstr(text_p);text_p = mynum2str(text_p,2);end	
%VB Macro
	%Selection.TypeText Text:="Test!"
	%in Matlab
	%set(word.Selection,'Text','test');
	%this also works
	%word.Selection.TypeText('This is a test');    
    if(enters_p(1))
        actx_word_p.Selection.TypeParagraph; %enter
    end
	%actx_word_p.Selection.Style = style_p;
    if(nargin == 6)%check to see if color_p is defined
        actx_word_p.Selection.Font.Color=color_p;     
    end
    
	actx_word_p.Selection.TypeText(text_p);
%    actx_word_p.Selection.Font.Color='wdColorAutomatic';%set back to default color
    for k=1:enters_p(2)    
        actx_word_p.Selection.TypeParagraph; %enter
    end
    actx_word_p.selection.Paragraphs.Alignment = sel0;
end
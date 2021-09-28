function Extra_code_S4(figure1, subplot1, subplot2)

xlim(subplot1,[0.0049141264204139 40625778.5455678]);
ylim(subplot1,[0.000660894157161765 0.0314005397084759]);

xlim(subplot2,[0.0287860005734426 813274.595671752]);
ylim(subplot2,[0.000215455548253829 0.0656210196600505]);


set(subplot2,'XMinorTick','on','XScale','log','YMinorTick','on','YScale',...
    'log');

% Create textboxs:---------------------------------------------------------
annotation(figure1,'textbox',...
    [0.403425578831311 0.424162257495587 0.0475115766262406 0.0458553791887129],...
    'VerticalAlignment','middle',...
    'String','b=n-q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.403425578831311 0.161375661375658 0.0475115766262407 0.0458553791887127],...
    'VerticalAlignment','middle',...
    'String','b<n-q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.143880926130099 0.650793650793649 0.156560088202867 0.0458553791887127],...
    'VerticalAlignment','middle',...
    'String','Increasing A_0',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

annotation(figure1,'textbox',...
    [0.403425578831311 0.75132098765432 0.0475115766262407 0.0458553791887126],...
    'VerticalAlignment','middle',...
    'String','b>n-q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(figure1,'textbox',...
    [0.751826901874309 0.583774250440917 0.143432194046308 0.0705467372134018],...
    'VerticalAlignment','middle',...
    'String','Increasing q',...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create arrows:-----------------------------------------------------------
annotation(figure1,'arrow',[0.305402425578831 0.305402425578831],...
    [0.619811618165782 0.756614087301586]);

annotation(figure1,'arrow',[0.751929437706725 0.751929437706725],...
    [0.310405643738977 0.714285714285714]);

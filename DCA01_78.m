clc;
clear all;
close all;
%data = 'aabaacabcabcb';
data = '00000011010100000110101';
%data = 'ABBABBABBBAABABA';
disp("Message :");
disp(data);
n=length(data);
parse=data(1);
i=2;
temp='';
encode='';
encode=strcat(encode,'(');
encode=strcat(encode,'0');
encode=strcat(encode,',');
encode=strcat(encode,parse);
encode=strcat(encode,')');
temp=strcat(temp,data(i));
prev=0;
while i<=n
    flag=0;
    C=strsplit(parse,',');
    for j=1:length(C)
        if strcmp(C(j),temp)
            flag=1;
            prev=j;
            break;
        end
    end
    if flag==1
        i=i+1;
        if i<=n
            temp=strcat(temp,data(i));
        else
            parse=strcat(parse,',');
            parse=strcat(parse,temp);
            encode=strcat(encode,',');
            encode=strcat(encode,'(');
            encode=strcat(encode,int2str(prev));
            encode=strcat(encode,',');
            encode=strcat(encode,'$');
            encode=strcat(encode,')');
        end
    else
        parse=strcat(parse,',');
        parse=strcat(parse,temp);
        encode=strcat(encode,',');
        encode=strcat(encode,'(');
        encode=strcat(encode,int2str(prev));
        encode=strcat(encode,',');
        encode=strcat(encode,temp(length(temp)));
        encode=strcat(encode,')');
        prev=0;
        i=i+1;
        if i<=n
            temp=data(i);
        end
    end
end
disp('Parsed message :');
disp(parse);
disp('Encoded message :');
disp(encode);
encode = strip(encode,'(');
encode = strip(encode,'right',')');
C4 = strsplit(encode,'),(');
track='';
decode='';
for i=1:length(C4)
    temp=C4(i);
    C5=strsplit(char(temp),',');
    if strcmp(C5(2),'$')
        decode=strcat(decode,'$');
    elseif strcmp(C5(1),'0')
        decode=strcat(decode,char(C5(2)));
        if i~=1
            track=strcat(track,',');
        end
        track=strcat(track,char(C5(2)));
    else 
        index=round(str2double(C5(1)));
        r1=strsplit(track,',');
        r2='';
        r2=strcat(r2,char(r1(index)));
        r2=strcat(r2,char(C5(2)));
        decode=strcat(decode,r2);
        track=strcat(track,',');
        track=strcat(track,r2);
    end
end
disp("Decoded message :");
disp(decode);
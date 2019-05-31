clc;
clear all;
close all;
%data = 'aacaacabcabaaac';
%data = '00000011010100000110101';
data = 'ABBABBABBBAABABA';
disp("Message :");
disp(data);
i = 1;
n = length(data);
parse = '';
temp = '';
flag=0;
prev1=1;
flag_val='';
pointer='';
length_val='';
w=4;
while i<=n
    prev=i;
    j=i-w;
    if j<=0
        j=1;
    end
    k1=0;
    k2=0;
    flag=0;
    temp='';
    temp1='';
    while j<i && i<=n
        if data(i)==data(j)
            temp=strcat(temp,data(i));
            i=i+1;
            j=j+1;
            if i==n+1
                k1=i;
                k2=j;
                temp1=strcat(temp1,temp);
            end
            flag=1;
        elseif data(i)~=data(j) && flag==1
            flag=0;
            if length(temp1)<length(temp)
                temp1='';
                temp1=strcat(temp1,temp);
                temp='';
                k1=i;
                k2=j;
                prev1=i;
            end
            if j<prev
                i=prev;
            else
                i=prev1;
                break;
            end
        else
            j=j+1;
        end
    end
    if length(temp1)==0 && i<=n
        temp1='';
        temp1=strcat(temp1,data(i));
        flag_val=strcat(flag_val,'0');
        flag_val=strcat(flag_val,',');
        i=i+1;
    else
        flag_val=strcat(flag_val,'1');
        flag_val=strcat(flag_val,',');
        pointer=strcat(pointer,int2str(k1-k2));
        pointer=strcat(pointer,',');
        length_val=strcat(length_val,int2str(length(temp1)));
        length_val=strcat(length_val,',');
    end
    if i==j
        i=prev1;
    end
    parse=strcat(parse,temp1);
    if i<=n
        parse=strcat(parse,',');
    end
end
disp('Parsed message :');
disp(parse);
%disp(flag_val);
%disp(pointer);
%disp(length_val);
C1 = strsplit(parse,',');
C2 = strsplit(pointer,',');
C3 = strsplit(length_val,',');
j=1;
k=1;
temp2="";
for i=1:length(C1)
    if flag_val(k)=='0'
        temp2=strcat(temp2,'(');
        temp2=strcat(temp2,flag_val(k));
        temp2=strcat(temp2,',');
        temp2=strcat(temp2,C1(i));
        temp2=strcat(temp2,')');
        k=k+2;
    else
        temp2=strcat(temp2,'(');
        temp2=strcat(temp2,flag_val(k));
        temp2=strcat(temp2,',');
        temp2=strcat(temp2,C2(j));
        temp2=strcat(temp2,',');
        temp2=strcat(temp2,C3(j));
        temp2=strcat(temp2,')');
        j=j+1;
        k=k+2;
    end
    if i<length(C1)
        temp2=strcat(temp2,',');
    end
end
disp('Encoded message :');
disp(temp2);
temp2 = strip(temp2,'(');
temp2 = strip(temp2,'right',')');
C4 = strsplit(temp2,'),(');
decode='';
for i=1:length(C4)
    temp=C4(i);
    C5=strsplit(temp,',');
    if strcmp(C5(1),'0')
        decode=strcat(decode,char(C5(2)));
    else
        point = str2double(C5(2));
        l = str2double(C5(3));
        for j=(length(decode)-round(point)+1):(length(decode)-round(point)+l)
            decode=strcat(decode,decode(j));
        end
    end
end
disp("Decoded message :");
disp(decode);
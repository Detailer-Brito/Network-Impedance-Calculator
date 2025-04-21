import numpy as np
import math
from sympy import *
import cv2
print(cv2.__version__)
from PIL import Image, ImageTk
import os

IMAGE_SIZE = (250, 250)  # Resize images to fit the frame


def insertImage(Imp_C,Imp_T,obj1,val2,config):
    print(config)
    image=cv2.imread('Images/'+config+'.png')
    array_coord=[(126,34),(184,110),(122,104),(227,36)]
    if(config=='LPLS' or config=='CPLS' or config=='LPCS' or config=='CPCS'):
        coord1=array_coord[2]
        coord2=array_coord[3]
    else:
        coord1=array_coord[0]
        coord2=array_coord[1]

    image=cv2.putText(image, "{:.2E}".format(val2[0]),coord1,cv2.FONT_HERSHEY_DUPLEX,0.4,(0,0,0))
    image=cv2.putText(image, "{:.2E}".format(val2[1]),coord2,cv2.FONT_HERSHEY_DUPLEX,0.4,(0,0,0))

    cv2.imwrite("image_dump.png",image)
    img_path = os.path.join("image_dump.png")
    img = Image.open(img_path)
    img.thumbnail(IMAGE_SIZE)
    image_put = ImageTk.PhotoImage(img)
    os.remove("image_dump.png")


    return image_put

    #obj1.setPixmap(QtGui.QPixmap.fromImage(image))
    #Saber width e length do objeto

def algoritmo(sreal,simag,creal,cimag,zcara,freq,obj):
    
    images_display=[]

    frequencia=freq
    print(cimag)
    Impedance_Linha=zcara
    Impedance_Carga=np.complex(creal,cimag)
    Impedance_T=np.complex(sreal,simag)
    Impedance_Norm_Carga=Impedance_Carga/Impedance_Linha
    Impedance_Norm_T=Impedance_T/Impedance_Linha
    condutancia_carga=np.real(1/Impedance_Norm_Carga)
    impedancia_carga=np.real(Impedance_Norm_Carga)

    R_Carga=np.real(Impedance_Norm_Carga)
    C_Carga=np.real(1/Impedance_Norm_Carga)
    R_T=np.real(Impedance_Norm_T)
    C_T=np.real(1/Impedance_Norm_T)
    Ks_T=(Impedance_Norm_T-1)/(1+Impedance_Norm_T)
    angle_inicio=np.angle(Ks_T)
    angulos2=np.linspace(angle_inicio,math.pi,10000)
    angulos3=np.linspace(-math.pi,angle_inicio,10000)
    angulos4=np.concatenate((angulos2,angulos3))
    pontos_rt=(1/(1+R_T))*np.cos(angulos4)+(R_T/(1+R_T)) + 1j*(1/(1+R_T))*np.sin(angulos4)
    pontos_ct=(1/(1+C_T))*np.cos(angulos4)+(-C_T/(1+C_T)) + 1j* (1/(1+C_T))*np.sin(angulos4)
    pontos_marcado=0
    pontos_marcado_2=0
    for i in range(0,np.size(pontos_rt)):

            condutancia=np.real((1-pontos_rt[i])/(1+pontos_rt[i]))
            if(condutancia>condutancia_carga-0.002 and condutancia <condutancia_carga+0.002):
                pontos_marcado=pontos_rt[i];
            
            
            impedancia=np.real((1+pontos_ct[i])/(1-pontos_ct[i]))
            if(impedancia>impedancia_carga-0.002 and impedancia <impedancia_carga+0.002):
                pontos_marcado_2=pontos_ct[i];

    config1=["","","",""]
    pontos_marcado_condut=np.conj(pontos_marcado)
    pontos_array=[pontos_marcado_condut,pontos_marcado]
    pontos_array_2=[pontos_marcado_2,np.conj(pontos_marcado_2)]
    valores=np.zeros((4,2))
    if(pontos_marcado!=0):
        for i in range(0,np.size(pontos_array)):
            aux=pontos_array[i];
            imped_dif=np.imag((1+aux)/(1-aux))-np.imag(Impedance_Norm_T)
            if(imped_dif>0):
                #Tem-se um indutor em série
                L=(imped_dif*Impedance_Linha)/(2*math.pi*frequencia);
                valores[i,0]=L;
                config1[i]="LS";
            else:
                #Tem-se um condensador em série
                C=np.abs((1/(2*math.pi*frequencia*imped_dif*Impedance_Linha)));
                valores[i,0]=C;
                config1[i]="CS";
            

            admit_dif=np.imag((1-aux)/(1+aux))-np.imag(1/np.conj(Impedance_Norm_Carga));
            if(admit_dif>0):
                #Tem-se um indutor em paralelo
                L=np.abs(Impedance_Linha/(2*math.pi*frequencia*admit_dif));
                valores[i,1]=L;
                config1[i]=config1[i]+"LP";
            else:
                #Tem-se um condensador em paralelo
                C=np.abs(admit_dif/(2*math.pi*frequencia*Impedance_Linha));
                valores[i,1]=C;
                config1[i]=config1[i]+"CP";

            image_aux=insertImage(Impedance_Carga,Impedance_T,obj,valores[i,:],config1[i])
            print(image_aux)
            images_display.append(image_aux)
    else:
        i=-1

    if(pontos_marcado_2!=0):
        for j in range(0,np.size(pontos_array_2)):
            aux=pontos_array_2[j]
            
            admit_dif=np.imag((1-aux)/(1+aux))-np.imag(1/Impedance_Norm_T);
            if(admit_dif<0):
                #Tem-se um indutor em paralelo
                L=np.abs(Impedance_Linha/(2*math.pi*frequencia*admit_dif));
                valores[i+j+1,0]=L
                config1[i+j+1]="LP"
            else:
                #Tem-se um condensador em paralelo
                C=admit_dif/(2*math.pi*frequencia*Impedance_Linha);
                valores[i+j+1,0]=C;
                config1[i+j+1]="CP";
            
            
            imped_dif=np.imag((1+aux)/(1-aux))-np.imag(np.conj(Impedance_Norm_Carga));
            if(imped_dif<0):
                #Tem-se um indutor em série
                L=np.abs((imped_dif*Impedance_Linha)/(2*math.pi*frequencia));
                valores[i+j+1,1]=L;
                config1[i+j+1]=config1[i+j+1]+"LS";
            else:
                #Tem-se um condensador em série
                C=np.abs((1/(2*math.pi*frequencia*imped_dif*Impedance_Linha)));
                valores[i+j+1,1]=C;
                config1[i+j+1]=config1[i+j+1]+"CS";

            image_aux=insertImage(Impedance_Carga,Impedance_T,obj,valores[i+j+1,:],config1[i+j+1])
            print(image_aux)
            images_display.append(image_aux)
    
    return images_display




def insertImage_stubs(Imp_C,Imp_T,obj1,val2,config):
    print(config)
    image=cv2.imread('Images/'+config+'.png')
    array_coord=[(188,50),(131,127),(82,79)]
    if(config=='SB_CC_P' or config=='SB_CA_P' ):
        coord1=array_coord[0]
        coord2=array_coord[1]
    else:
        coord1=array_coord[0]
        coord2=array_coord[2]

    aux="{:.3E}".format(val2[0])#str(val2[0])
    image=cv2.putText(image, "l="+aux,coord1,cv2.FONT_HERSHEY_DUPLEX,0.4,(0,0,0))
    aux="{:.3E}".format(val2[1])#str(val2[1])
    image=cv2.putText(image, "l="+aux,coord2,cv2.FONT_HERSHEY_DUPLEX,0.4,(0,0,0))

    cv2.imwrite("image_dump.png",image)
    img_path = os.path.join(IMAGE_DIR, "image_dump.png")
    img = Image.open(img_path)
    img.thumbnail(IMAGE_SIZE)
    image_put = ImageTk.PhotoImage(img)
    os.remove("image_dump.png")


    return image_put


def algoritmo_stubs(sreal,simag,creal,cimag,zcara,freq,obj):

    images_display=[]

    frequencia=freq
    Impedance_Linha=zcara
    Impedance_Carga=np.complex(creal,cimag)
    Impedance_T=np.complex(sreal,simag)

    Impedance_Norm_Carga=Impedance_Carga/Impedance_Linha
    condutancia_carga=np.real(1/Impedance_Norm_Carga)
    impedancia_carga=np.real(Impedance_Norm_Carga)
    Impedance_Norm_T=Impedance_T/Impedance_Linha
    Ks_carga=(Impedance_Norm_Carga-1)/(1+Impedance_Norm_Carga)
    Ks_T=(Impedance_Norm_T-1)/(1+Impedance_Norm_T)
    angle_inicio=np.angle(Ks_T)
   


    R_Carga=np.real(Impedance_Norm_Carga)
    C_Carga=np.real(1/Impedance_Norm_Carga)
    R_T=np.real(Impedance_Norm_T)
    C_T=np.real(1/Impedance_Norm_T)
    

    angle_inicio=np.angle(Ks_carga)

    if(angle_inicio>0):
        angulos2=np.linspace(angle_inicio,0,10000)
        angulos3=np.linspace(0,-2*math.pi,10000)
        angulos4=np.linspace(angle_inicio-math.pi,-math.pi,10000)
        angulos5=np.linspace(-math.pi,-2*math.pi,10000)
    else:
        angulos2=np.linspace(angle_inicio,-math.pi,10000)
        angulos3=np.linspace(-math.pi,-(2*math.pi-angle_inicio),10000)
        angulos4=np.linspace(angle_inicio+math.pi,0,10000)
        angulos5=np.linspace(0,-(2*math.pi-angle_inicio+math.pi),10000)
    


    #Circunferencia com raio Ks
    pontos=np.abs(Ks_carga)*np.cos(np.concatenate((angulos2,angulos3))) + 1j*np.abs(Ks_carga)*np.sin(np.concatenate((angulos2,angulos3)));
    pontos_a=np.abs(Ks_carga)*np.cos(np.concatenate((angulos4,angulos5))) + 1j* np.abs(Ks_carga)*np.sin(np.concatenate((angulos4,angulos5)));
    

    swr=(1+np.abs(Ks_carga))/(1-np.abs(Ks_carga))
    cho=[]
    chose=[]
    #Condicao de adaptacao
    if(np.real(Impedance_Norm_T)>np.real(swr) or np.real(Impedance_Norm_T)<np.real(1/swr)):
        print("Nao e possivel fazer adaptacao por stubs")
    else:
        for i in range(0,np.size(pontos)): 
            
            res=np.real((1+pontos[i])/(1-pontos[i]))
            
            if(res>R_T-0.002 and res<R_T+0.002):
                cho=np.append(cho,pontos[i])#[cho pontos(i)]
                
          
            
            admit=np.real((1+pontos_a[i])/(1-pontos_a[i]))
            if(admit>C_T-0.002 and admit<C_T+0.002):
                chose=np.append(chose,pontos_a[i])#[chose pontos_a(i)]
                
      
        
        config1=["","","",""]
        valores=np.zeros((4,2))
        ponto_chosen=cho[1]
        ponto_admit=chose[1]
        pontos_array=[ponto_chosen,np.conj(ponto_chosen)]
        pontos_array_2=[ponto_admit,np.conj(ponto_admit)]
        #Calculo do stub em série

        if(ponto_chosen!=0):
            for i in range(0,np.size(pontos_array)):
                #Calculo da distancia de linha
                aux=pontos_array[i];
                ang_linha=np.angle(aux)-np.angle(Ks_carga);
                if(ang_linha>0):
                    ang_linha=2*math.pi-ang_linha;
                else:
                    ang_linha=abs(ang_linha);
                

                factor_linha=4*math.pi/ang_linha;
                distancia=(3*10**8)/(frequencia*factor_linha);
                valores[i,0]=distancia;
                im_dif=np.imag(np.conj(Impedance_Norm_T))-np.imag((1+aux)/(1-aux))
                
                if(im_dif>0):
                    #Stub CC
                    factor=2*math.pi/np.arctan(im_dif);
                    l=((3*10**8)/(frequencia))/factor;
                    config1[i]="SB_CC_S";
                else:
                    #Stub CA
                    factor=2*math.pi/acot(im_dif);
                    l=np.abs(((3*10**8)/(frequencia))/factor);
                    config1[i]="SB_CA_S";

                valores[i,1]=l;
                images_display.append(insertImage_stubs(Impedance_Carga,Impedance_T,obj,valores[i,:],config1[i]))

        #Calculo do stub em paralelo
        
        if(ponto_admit!=0):
            for j in range(0,np.size(pontos_array_2)):

                aux=pontos_array_2[j];
            
                ang_linha_admit=np.angle(aux)-(np.angle(Ks_carga)+math.pi);
                if(ang_linha_admit>0):
                    ang_linha_admit=2*math.pi-ang_linha_admit;
                else:
                    ang_linha_admit=abs(ang_linha_admit);
                
                factor_linha=4*math.pi/ang_linha_admit;
                distancia_admit=(3*10**8)/(frequencia)*1/factor_linha;
                valores[i+j+1,0]=distancia_admit
                admit_dife=np.imag(1/np.conj(Impedance_Norm_T))-np.imag((1+aux)/(1-aux));
                if(admit_dife<0):
                    #Stub CC
                    factor_admit=2*math.pi/np.arctan(-1/admit_dife)
                    l_admit=((3*10**8)/(frequencia))/factor_admit;
                    config1[i+j+1]="SB_CC_P"
                    
                else:
                    #Stub CA
                    factor_admit=2*math.pi/acot(1/admit_dife)
                    l_admit=((3*10**8)/(frequencia))/factor_admit;
                    config1[i+j+1]="SB_CA_P"

                valores[i+j+1,1]=l_admit
                images_display.append(insertImage_stubs(Impedance_Carga,Impedance_T,obj,valores[i+j+1,:],config1[i+j+1]))
        
        return images_display
        print(valores)
        print(config1)
    
Program Extensao;
(*------------------------------------------------------------------------
  Programa SEQEx - Simulacao de Equilibrios por Extensao

  Desenvolvido por:  Prof.Dr.Claudimir Lucio do Lago
                     Robinson Naohide Mioshi

  Inicio do projeto: 30/11/94
-------------------------------------------------------------------------*)
uses dos,crt,lu1,geral,graph;

const
  MAXESP = 40;                     { Numero maximo de especies }
  MAXEQ = 30;                      { Numero maximo de equilibrios }
  MAXFASE = 4;                     { Numero maximo de fases }

type
  TipoEspecie = Record
    nmol : Float;                  { Numero de mols }
    nome : String;                 { Especie }
    C    : Float;                  { Concentracao da especie mol/L }
    carga: Integer;                { carga da especie }
 diametro: Float;                  { diametro efetivo da especie }
    fase : Integer;                { fase que pertence a especie }
  end;
  TipoListaEspecie = Record
    nesp : Byte;                              { numero total de especie }
    esp  : Array[1..MAXESP] of TipoEspecie;   { numero da especie }
  end;
  PEspecie = ^TipoListaEspecie;
  TipoSistema = Record
    K,                      { constante de equilibrio termodinamico }
    K_1,
    fm : Float;
  end;
  TipoFase = (pura,mistura);         { tipo de fase }
  TipoListaSistema = Record
    neq : Integer;                              { Numero de equilibrio }
    equ : Array[1..MAXEQ] of TipoSistema;
   coef,                                        { coefiente estequiometrico }
   expo: Array[1..MAXESP,1..MAXEQ] of ShortInt; { expoente do Quociente     }
    vol : Array[0..MAXFASE] of Float;           { volume da fase }
    tf  : Array[0..MAXFASE] of TipoFase;        { numero da fase }
  end;
  PSistema = ^TipoListaSistema;

Procedure LeReacao(   nomearq : String;
                  var sistema : PSistema;
                  var especie : PEspecie );
(**************************************************************************
Objetivo : Ler cada equilibrio do arquivo de reacao (d_XXXX.dat)
Entrada  : nomearq - Nome do arquivo
           sistema - Dados da especie no equilibrio
           especie - Dados apena da especie
Saida    : sistema
           especie
Erro     :
**************************************************************************)
var arq : Text;
      s : String;     { Contem a reacao e a constante de equilibrio }
    pont: Byte;       { Indica a posicao da leitura numa linha }
    codigo,
    i,j : Integer;
    branco : Boolean;

Procedure LeEspecie(    aux_coef     :Integer;
                    var numeroesp    :Integer);
(**************************************************************************
Objetivo:  Le cada especie das reacoes
           nome da especie, carga da especie
Entrada :  aux_coef  - carga da especie
Saida   :  numeroesp - numero total de especie
Erro    :
**************************************************************************)
var                             c    : String;
            tamanho, codigo, i,carga : integer;
                          inicio,fim : Byte;
                           aux_carga : Float;
Begin
  inicio := pont;
  ProcuraCar(s,pont,' ',codigo);
  tamanho:=pont-inicio;
  If (codigo=0) and (tamanho>=3) then
  begin
    c:=Copy(s,inicio,tamanho);
    With especie^  do
    begin
      i:=1;
      While(i<=nesp) and (c<>especie^.esp[i].nome) do Inc(i);
      If i>nesp then
        begin
          Inc(nesp);
          esp[nesp].nome:=c;            { nome da especie }
          { Le carga da especie }
          inicio:=1;
          ProcuraConjCar(c,inicio,'+-',codigo);
          If codigo=0 then
            begin
              fim:=1;
              ProcuraCar(c,fim,'}',codigo);
              If fim=(inicio+1) then
                begin
                  Val(c[inicio]+'1',aux_carga,codigo);
                  carga:=Trunc(aux_carga);
                end else LeInt(c,inicio,carga,codigo);
            end else
          If codigo=2 then carga:=0;
          numeroesp:= i;
          esp[nesp].carga:=carga;             { carga da especie }
          { Le fase da especie }
          ProcuraCar(c,inicio,'^',codigo);
          If codigo = 0 then
            begin
              Inc(inicio);
              LeInt(c,inicio,esp[nesp].fase,codigo);
            end
          else esp[nesp].fase:=0;
        end
      else
      numeroesp:=i;
    end;
  end;

End;

Function LeComponente:Boolean;
(**************************************************************************
Objetivo: Leitura do coeficiente estequiometrico, constante de equilibrio
          a partir de cada reacao.
Entrada : s - equacao quimica
Saida   : sistema.coef[]
          sistema.equ[].K
Erro    :
**************************************************************************)
var   c                       : String;
    inicio,tamanho,codigo,i,
    aux_coef,numeroesp,fator  : Integer;
    fim                       : Boolean;
Begin
  fator:=-1;
  fim:= False;
  Repeat
    IgnoraCar(s,pont,' ',codigo);
    Case codigo of
    0 : Case s[pont] of
        '0'..'9','{':
                      begin
                        If s[pont]='{' then  aux_coef:=1 else
                          LeInt(s,pont,aux_coef,codigo);
                        LeEspecie(aux_coef,numeroesp);
                        { coeficiente estequiometrico da especie }
                        sistema^.coef[numeroesp,sistema^.neq]:=(fator*aux_coef);
                      end;
                '+' :Inc(pont);
                 '=':begin
                       fator:=1;

                       Inc(pont);
                     end;
            'K','k' :begin
                        Inc(pont);
                        LeReal(s,pont,sistema^.equ[sistema^.neq].K,codigo);
                        { Inversao do K }
                        sistema^.equ[sistema^.neq].K_1:=1/sistema^.equ[sistema^.neq].K;
                        sistema^.equ[sistema^.neq].fm:=sistema^.equ[sistema^.neq].K_1;
                        fim:=True;
                     end;
        end;
    end;
  Until fim;
End;
Begin   { LeReacao }

  {Zera os diametros }
  For i:= 1 to MAXESP do with especie^.esp[i] do
  begin
    diametro:=0.0;
    nmol:=0.0;
    C:=0.0;
  end;

  sistema^.neq:=1;
  especie^.nesp:=0;
  With sistema^ do
  { Zerar a matriz dos coeficientes }
  For i:= 1 to MAXEQ do
    For j:= 1 to MAXESP do
      sistema^.coef[j,i]:=0;
  branco:=true;;
  If AbrirTexto(arq,nomearq)  then
    While ((not Eof(arq)) and (branco=true))  do
      begin
        pont:=1;
        Readln(arq,s);
        If s<>'' then
          begin
            If (LeComponente)  then Inc(sistema^.neq);
          end
        else branco:=false;
      end;
  { Retirar }
  sistema^.neq:=sistema^.neq-1;
  Close(arq);
End;

Procedure LeConcentr(nomearq     : String;
                     var especie : PEspecie);
(**************************************************************************
Objetivo:  Le o numero de mols das especies do arquivo (c_XXXX.dat)
Entrada :  nomearq - Nome do arquivo
           especie - Dados das especies presentes no sistema.
Saida   :  Dados referentes ao numero de mols da especie, carga da especie e
           numero total de especies ( presente ou nao nas reacoes )
Erro    :
**************************************************************************)

var  i   : Integer;
     arq : text;
     s   : String;
   branco: Boolean;

Function LeEspecie(var especie : TipoListaEspecie  ): Boolean;
(**************************************************************************
Objetivo:  Le o numero de mols da especie, e ser for nova, incluir na lista
           ( carga, numero de mols).
Entrada :  especie : TipoListaEspeci - Dados das especies.
Saida   :  LeEspecie - true  ( Especie presente ).
                       false ( Especie nao presente ).
Erro    :
**************************************************************************)
var             fim,inicio : Byte;
    i,code,tamanho, codigo,carga : Integer;
                     pont  : Byte;
                         c : String;
        aux_nmol,aux_carga : Float;
Begin

  LeEspecie := true;
  pont:= 1;
  ProcuraCar(s,pont,'{',code);
  If code = 0  then
    begin
      inicio := pont;
      ProcuraCar(s,pont,' ',code);
      If code = 0 then with especie do
        begin
          tamanho := pont - inicio;
          c := Copy(s,inicio,tamanho);
          i:= 1;
          While ( c <> esp[i].nome ) and (i <= nesp)  do Inc(i);
          If i <= nesp  then
          begin
            LeReal(s,pont,aux_nmol,code);
            If code = 0 then esp[i].nmol := esp[i].nmol + aux_nmol;
            ProcuraCar(s,pont,' ',code);
            If (code=2) or (pont=length(s)) then esp[i].diametro:=0.0
            else
              begin
                LeReal(s,pont,esp[i].diametro,code);
                  If code=2 then  esp[i].diametro := 0.0;
              end;
          end
          else
            begin
              nesp := i;
              { Concentracao da especie }
              esp[i].nome := c;

              { Le carga da especie }
              inicio:=1;
              ProcuraConjCar(c,inicio,'+-',codigo);
              If codigo=0 then
                begin
                  fim:=1;
                  ProcuraCar(c,fim,'}',codigo);
                  If fim=(inicio+1) then
                    begin
                      Val(c[inicio]+'1',aux_carga,codigo);
                      carga:=Trunc(aux_carga);
                    end else LeInt(c,inicio,carga,codigo);
                end else
              If codigo=2 then carga:=0;
              { carga da especie }
              esp[nesp].carga:=carga;

              { Le fase da especie }
              ProcuraCar(s,inicio,'^',codigo);
              If codigo = 0 then
                begin
                  Inc(inicio);
                  LeInt(c,inicio,esp[nesp].fase,codigo);
                end
              else esp[nesp].fase:=0;

              ProcuraCar(s,pont,'=',code);
              ProcuraConjCar(s,pont,'-+0123456789Ee.',code);
              LeReal(s,pont,aux_nmol,code);
              If code = 0 then esp[i].nmol := aux_nmol;

              { Le o diametro do ion}
              esp[nesp].diametro:=0.0;
              LeReal(s,pont,esp[nesp].diametro,code);
              If code=2 then  esp[nesp].diametro := 0.0;
            end;
        end else LeEspecie:=false;
    end;
End;

Begin {LeConcentracao}
  { zera as concentracoes de todas as especies }
  With especie^ do
    For i:= 1 to nesp do with esp[i] do
      nmol := 0.0;
  { Faz a leitura da concentracao de todas as especies }
  If AbrirTexto(arq,nomearq) then
    begin
      branco:=false;
      While ((not Eof(arq)=true) and (branco=false)) do
        begin
          Readln(arq,s);
          If s<>'' then
          begin
            If not(LeEspecie(especie^)) then Writeln('Erro read file number mols');
          end else branco:=true;
        end;
    end;
  Close(arq);
End;

Procedure LeVolume (nomearq       : String;
                    var sistema   : PSistema;
                    var especie   : PEspecie);
(**************************************************************************
Objetivo: Le o volume de cada fase em mL
          MX :  Fase mistura X.
          PY :  Fase pura    Y.
Entrada : sistema : PSistema - Dados do equacao
          especie : PEspecie - Dados das especies.
Saida   : volume de cada fase
          fase pura = volume igual a 1.
Erro    :
***************************************************************************)

var
             arq : text;
               s : string;
  i,j,code, fase : integer;
            pont : Byte;
           tfase : TipoFase;
          branco : Boolean;
Begin
  { Zerar os volumes de todas fases }
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    sistema^.vol[fase]:=0.0;


  branco:=false;
  If AbrirTexto(arq,nomearq) then
  begin
    While ((not Eof(arq)) and (branco=false)) do
    begin
      Readln(arq,s);
      If s<>''  then
        begin
           pont:=1;
           ProcuraConjCar(s,pont,'MPmp',code);
             If code = 0 then
               begin
                 If UpCase(s[pont]) = 'M' then tfase := mistura
                   else tfase := pura;
                 ProcuraConjCar(s,pont,'0123456789',code);
                 If code = 0 then
                   begin
                     LeInt(s,pont,fase,code);
                     sistema^.tf[fase] := tfase;
                     If tfase = mistura then
                       begin
                         ProcuraConjCar(s,pont,'0123456789',code);
                         If code = 0 then LeReal(s,pont,sistema^.vol[fase],code);
                       end
                       else
                         sistema^.vol[fase] := 1.0;
                   end;
               end;
        end
      else branco:=true;
    end;
  end;
  Close(arq);

  { expo }
  sistema^.expo:=sistema^.coef;
  For i:= 1 to sistema^.neq do with sistema^ do
    For j:= 1 to especie^.nesp do
    begin
      If tf[especie^.esp[j].fase]=pura then expo[j,i]:=0;
     end;
End;


Procedure Raio_Efetivo( especie :PEspecie);
var         arq : Text;
    inicio, fim : Byte;
    s,aux_s     : String;
    r           : Float;
    erro,i,pont : Integer;
Begin
  Assign(arq,'c:\tp\dadosn\raio.dat');
  Reset(arq);
  While (not Eof(arq)) do
    begin
      Read(arq,r);
      inicio:=1;
      fim:=1;
      pont:=1;
      Readln(arq,aux_s);
      While pont<Length(aux_s) do
        begin
          { Procura especie }
          ProcuraCar(aux_s,inicio,'{',erro);
          ProcuraCar(aux_s,fim,'}',erro);
          s:=Copy(aux_s,inicio,((fim+1)-inicio));
          For i:= 1 to especie^.nesp do with especie^.esp[i] do
            If s=nome then diametro:=r;
          inicio:=fim;
          Inc(fim);
          pont:=fim;
        end;
    end;
end;



Procedure Espalha(var especie : PEspecie;
                  var sistema : PSistema);
(**************************************************************************
Objetivo: Distribuir os numeros de mols das especies iniciais.
          Utiliza o processo da bisseccao para distribuir simetricamente
          os numeros de mols do reagente e produto.
Entrada : especie - Dados das especies.
          sistema - Dados da equacao de reacao.
Saida   : Numero de mols ditribuidos para reagentes e produtos.
Erro    :
**************************************************************************)


type equilibrio = Array[1..MAXEQ] of Float;

var j,i                   : Integer;
    psi_max,
    psi_minimo            : equilibrio;
    f : equilibrio;
    psi_H,psi_L,psi_zero  : Float;
                     cont : Boolean;

Function quociente_reacao( n:integer): Float;
var aux_numerador,aux_denominador : Float;
    i,j                           :integer;
Begin
    For j := 1 to especie^.nesp do
      begin
        aux_denominador:=1.0;aux_numerador:=1.0;
        If sistema^.expo[j,n]<0 then with especie^.esp[j] do
           aux_denominador:=aux_denominador*Elevar((nmol/sistema^.vol[fase]),sistema^.coef[j,n])
        else
        If sistema^.expo[j,n]>0 then with especie^.esp[j] do
            aux_numerador  :=aux_numerador *Elevar((nmol/sistema^.vol[fase]),sistema^.coef[j,n]);
      end;
    quociente_reacao:=(aux_numerador/aux_denominador);
End;

Procedure Intervalo( i:integer; var psi_max,psi_minimo : equilibrio );
(**************************************************************************
Objetivo:  Determina para cada equilibrio o psi_maximo e psi_minimo. Utilizando
           um valor meio de psi para dar a extensao de reacao inicial, assim
           distribuindo simetricamente os numeros de mols.
Entrada :  numero de mols das especies.
           coeficiente estequiometrico das esquecies em cada equilibrio.
Saida   :  Valor de psi_maximo e psi_minimo.
Erro    :
**************************************************************************)

var aux_psi_max,aux_psi_min : Float;
                          j : integer;
Begin
  j:= 1;
  psi_max[i]:=100;
  psi_minimo[i]:=-100;
  { Minimo valor de psi_maximo }
  For j:= 1 to especie^.nesp do with sistema^ do
  begin
    If (coef[j,i]<>0) then with especie^ do
    begin
      If (coef[j,i]<0.0) then
      begin
        aux_psi_max:= -((esp[j].nmol)/(coef[j,i]));
        If (aux_psi_max < psi_max[i]) then psi_max[i]:=aux_psi_max;
      end  else
      { Maximo valor de psi_minimo }
      begin
        aux_psi_min:= -((esp[j].nmol) /(coef[j,i]));
        If (aux_psi_min > psi_minimo[i]) then psi_minimo[i]:=aux_psi_min;
      end;
    end;
  end;
End;

Function NaoEspalha(n: Integer): Boolean;
(**************************************************************************
Objetivo:  Se tiver um dos produto ou reagente com o numero de mols
           diferente de zero, nao distribuir os numeros de mols.
Entrada :
Saida   :
Erro    :
**************************************************************************)

var  ok : Boolean;
      i : Integer;
Begin
  ok:=true;
  For i:= 1 to especie^.nesp do
    begin
      If sistema^.coef[i,n] <> 0 then
        If (especie^.esp[i].nmol=0.0) then ok:=false;
    end;
  NaoEspalha:=ok;
End;

Begin   { Espalha }
  Repeat
    i:=1;
    cont:= true;
    While ((sistema^.neq)>=i) do with sistema^ do
    begin
      If Not(NaoEspalha(i)) then
        begin
          { Obter os valores de psi_maximo e psi_minimo }
          Intervalo(i,psi_max,psi_minimo);
          If (psi_max[i]<>0.0) or (psi_minimo[i]<>0.0) then
          begin
            psi_H:=psi_max[i];
            psi_L:=psi_minimo[i];
            psi_zero:= (psi_L+psi_H)/2;
            f[i]:=(quociente_reacao(i) - (1/equ[i].fm));
            If f[i]<0.0 then psi_L:=psi_zero
            else psi_H:=psi_zero;
            For j:= 1 to especie^.nesp do with especie^.esp[j] do
              If (coef[j,i]<>0)  then
                 nmol:= nmol + (coef[j,i]*(psi_zero))
          end else cont:=false;
        end;
        Inc(i);
     end;
  Until cont=true;
End;


Procedure Equilibrar( var especie : PEspecie;
                      var sistema : PSistema;
                          sinal   : Boolean;
                          fi      : Boolean);
(**************************************************************************
Objetivo: Neste procedimento inclui a utilizacao da extensao de reacao com
          com o metodo da convergencia glogal junto com Newton-Raphson.
Entrada : especie - Dados da especie
          sistema - Dados da reacao
Saida   : especie - numero de mols das especies no equilibrio.
Erro    :
**************************************************************************)


type equilibrio = Array[1..MAXEQ] of Float;
     esp        = Array[1..MAXESP] of Float;

var
   quociente,
   qsi,
   aux_K,temp_K            : equilibrio;
   mi, temp_C, temp_nmol   : esp;
   i,j,m,passo,cont        : Integer;
   soma,r,valor, f_ionica  : Float;
   fx,p        : TNMatrix;
   f,temp_f,d  : TNVector;
   erro        : Byte;
   nie         : Word;
   arq1        : Text;
   b           : glindx;
   sair        : String;


Procedure Forca_Ionica_K ( var k_linha: equilibrio;
                           var forca_I: Float );
(**************************************************************************
Objetivo: Calcula o K em funcao da forca ionica do meio pela equacao de
          Debye-Huckel, atraves da constante termodinamica.
Entrada : K       - constante termodinamica.
          especie - calcular a forca ionica do meio.
Saida   : K_linha - constante condicional (I=constante).
Erro    :
**************************************************************************)

var                   i,j : Integer;
                   coef_a : equilibrio;
    denominador,numerador : Float;
                   ativ   : esp;
Begin
  { Calcular a forma ionica }
  soma:=0.0;
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    If carga<>0 then soma:=soma+(carga*carga)*(temp_nmol[i]/sistema^.vol[fase]);
  forca_I:=soma*0.5;
  { Calcula a coeficiente de atividade de cada especie ionica }
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    begin
      If carga<>0 then
      begin
        coef_a[i]:=(exp(2.3025851*(-0.5*(carga*carga)*sqrt(forca_I))/(1+ 0.33*diametro*sqrt(forca_I))));
      { Writeln(nome,' = ',coef_a[i]);}
      end
      else coef_a[i]:=1;
    end;
  {  Calcular o K linha }
  For i:= 1 to sistema^.neq do with sistema^ do
  begin
    numerador:=1.0;
    denominador:=1.0;
    For j:= 1 to especie^.nesp do with especie^ do
      If coef[j,i]>0 then numerador:=numerador*Elevar(coef_a[j],coef[j,i]) else
      If coef[j,i]<0 then denominador:=denominador*Elevar(coef_a[j],-coef[j,i]);
    { *** ATENCAO ***** }
    { O K_linha esta invertido, como o K termodinamico }
    K_linha[i]:= (numerador*sistema^.equ[i].fm)/(denominador);
  end;
End;

Procedure Quocient;
(**************************************************************************
Objetivo:   So calcula o quociente de reacao para cada equilibrio.
Entrada :
Saida   :
Erro    :
**************************************************************************)

var                   i,j : Integer;
    numerador,denominador : Float;

Begin
  For i:= 1 to sistema^.neq do with sistema^ do
  begin
    numerador:=1.0;
    denominador:=1.0;
    For j:= 1 to especie^.nesp do with especie^ do
      If expo[j,i]>0 then numerador:=numerador*Elevar((temp_nmol[j]/vol[esp[j].fase]),coef[j,i]) else
      If expo[j,i]<0 then denominador:=denominador*Elevar((temp_nmol[j]/vol[esp[j].fase]),-coef[j,i]);
      quociente[i]:=numerador/denominador;
  end;
End;


Function Positivo : Boolean;
(**************************************************************************
Objetivo:  Verifica-se temp_nmol e' positivo ou nao.

Entrada : temp_nmol - Numero de mols das especies.
Saida   : Positivo - true : Numero de mols >= 0
          Positivo - false: Numeor de mols < 0
Erro    :
**************************************************************************)

var ok : Boolean;
     i : Integer;
Begin
  ok:=True;
  For i:= 1 to especie^.nesp do
    If temp_nmol[i]<=0.0 then ok:=False;
  Positivo:=ok;
End;

Function Equilibrou : Boolean;
(**************************************************************************
Objetivo: Verifica-se o sistema alcancou ou nao o equilibrio.
          Atraves do parametro da funcao f. (f=Q/K -1)
          Onde o criterio de convergencia e' :
          f > 1E-13  Nao convergiu para solucao.
          f <= 1E-13 Convergiu para solucao.

          Outro criterio de convergencia do ponto de vista quimico e' quando
          a funcao f mantem constante por algumas iteracoes ( escolhemos
          10 iteracoes).
Entrada : f   - funcao f ( f=Q/K -1).
          nie - Numero de iteracoes onde f mantem constante.
Saida   : Equilibrou -  true : equilibrou o sistema.
          Equilibrou -  false: o sistema nao esta equilibro.
Erro    :
**************************************************************************)

var
    j     : Word;
    ok,ok1: Boolean;
    reacao: Array[1..20] of Integer;
    i     : Integer;

Begin
  ok:=True;
  Equilibrou:=false;
  i:=0;
  For j:= 1 to sistema^.neq do
    If (Abs(f[j])>1E-13) then
      begin
        ok:=false;
        Inc(i);
        reacao[i]:=j;
      end;
  {Writeln('i =',i,'     ');}
  If ok then Equilibrou:=ok else
    begin
      ok:=true;
      For j:=1 to i do
        begin
          If (Abs(temp_f[reacao[j]]-f[reacao[j]])>1E-6) then ok:=false;
          If ((quociente[reacao[j]])<=(1/temp_K[reacao[j]])) and ((quociente[reacao[j]])>1) then ok:=false;
          If ((quociente[reacao[j]])>=(1/temp_K[reacao[j]])) and ((quociente[reacao[j]])<1) then ok:=false;
        end;
      If ok then Inc(nie) else nie:=0;
      Equilibrou:=nie=20;
      {Writeln('nie =',nie);}
    end;
    temp_f:=f;
{   For j:= 1 to sistema^.neq do
      begin
        Writeln('f[',j,'] =',f[j]);
        Writeln('temp_f[',j,'] =',temp_f[j]);
      end;
     Writeln(nie);}


End;

Function Constante_K : Boolean;
(**************************************************************************
Objetivo: Verifica-se a constante condicional mantem constante. Este e'
          obtido dos valores da forcao ionica do meio, coeficiente de
          atividade (equacao Debye-Huckel) e constante termodinamica
          ( Forca ionica constante ).
Entrada : temp_K - constante condicional atual.
          aux_K  - constante condicional anterior.
Saida   : Constante_K - true : Equilbrou.
                      - false: Nao equilbriou.
Erro    :
**************************************************************************)

var
    j : Word;
    ok: Boolean;
Begin
  ok:=True;
  For j:= 1 to sistema^.neq do
    If Abs(temp_K[j]-aux_K[j])/(temp_K[j]+aux_K[j])<=1E-10 then
    Constante_K:=true else
    Constante_K:=false;
End;


(************************ PRINCIPAL *******************************)

Begin  { Equilibrar }
  { Guarda os valores de K original }
  For j:= 1 to sistema^.neq do
    begin
      qsi[j] := 0.0;
      temp_K[j]:=sistema^.equ[j].fm;
    end;
  For i:= 1 to sistema^.neq do temp_f[i]:=0.0;
  For j:= 1 to especie^.nesp do  temp_nmol[j]:=especie^.esp[j].nmol;

  { O numero de moles e igual ao volume da especie }
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    If (sistema^.tf[fase]=pura) then sistema^.vol[fase]:=temp_nmol[i];

  cont:=0;
  { Loop do efeito da forcao ionica }
  passo:=0;
  nie:=0;
  Repeat
    Inc(cont);

    { Guarda os valores de K temp anterior}
    For j:= 1 to sistema^.neq do  aux_K[j]:=temp_K[j];

    With sistema^ do
    Repeat


      If not(sinal) then
      begin
        Inc(passo);
        GotoXY(1,1);
        Writeln('iterations  =',passo);
      end;


      { Atualiza o numero de mols }
      For j:= 1 to especie^.nesp do especie^.esp[j].nmol:=temp_nmol[j];

      { Calcula o vetor funcao (constante do sistema de equacoes). }
      Quocient;
      For i:= 1 to sistema^.neq do with sistema^.equ[i] do
      begin
        f[i]:=(quociente[i]*temp_K[i])-1;
        If f[i]<0 then
          begin
            fm:=1/fm;
            temp_K[i]:=1/temp_K[i];
            aux_K[i]:=1/aux_K[i];
            quociente[i]:=1/quociente[i];
            f[i]:=(quociente[i]*temp_K[i])-1;
            For j:= 1 to especie^.nesp do
              begin
                coef[j,i]:=-coef[j,i];
                expo[j,i]:=-expo[j,i];
              end;
          end;
      end;

      { Monta o Jacobiano }
      For i:= 1 to neq do
        For m:= 1 to neq do
        begin
          soma:=0.0;
          For j:= 1 to especie^.nesp do
             soma:= soma + {expo[j,i]}coef[j,i]*coef[j,m]/temp_nmol[j];
          fx[i,m]:=quociente[i]*soma*temp_K[i];
        end;

      { Resolve o sistema de equacoes para determinar o passo (d) no NR. }
      Ludcmp(fx,neq,b,valor,erro);
      Lubksb(fx,neq,b,f,d);
      If erro<>0 then
        begin
          Writeln('Erro : Exit Program');
          Readln;
       end;
      { Estabelece o passo e atualiza os numeros de moles. }
      r:=2.0;
      Repeat
        r:=0.5*r;
        For j:= 1 to neq do qsi[j] := -r*d[j];
        { Atualizar os numeros de mols }
        For i:= 1 to especie^.nesp do with especie^.esp[i] do
          begin
            temp_nmol[i] := nmol;
            For j:= 1 to neq do
                temp_nmol[i] := temp_nmol[i] + coef[i,j] * qsi[j];
          end;
      Until Positivo;

      { O numero de moles e igual ao volume da especie }
      For i:= 1 to especie^.nesp do with especie^.esp[i] do
        If sistema^.tf[fase]=pura then
           If temp_nmol[i]<1E-500 then
             begin
               temp_nmol[i]:=1E-40;
               sistema^.vol[fase]:=temp_nmol[i];
             end else sistema^.vol[fase]:=temp_nmol[i];



     { Para ver }
(*
      For i:=1 to sistema^.neq do
      begin
        Writeln('f[',i,'] =',f[i]);
       {Writeln('temp_K[',i,'] =',temp_K[i],'  quociente[',i,'] =',quociente[i]);}
      end;
*)
      Writeln;
      Writeln('        Components     Concentration         Numbers Mols');
      For i:= 1 to especie^.nesp  do  with especie^.esp[i] do
      Writeln(especie^.esp[i].nome:15,' =',temp_nmol[i]/sistema^.vol[fase],especie^.esp[i].nmol{,'  ',sistema^.vol[fase]});


    Until Equilibrou;

    If not(FI) then
      begin
        { Calcula a forma ionica com as concentracoes do equilibrio }
        Forca_Ionica_K (temp_K,f_ionica);

{       Writeln('forca ionica =',f_ionica);}

(*
        { ver }
        For i:= 1 to sistema^.neq do
          begin
            Writeln('K=',sistema^.equ[i].fm,' K linha= ',(temp_K[i]));
            Writeln('K_aux =',aux_K[i]);
          end;
*)
      end;

  Until Constante_K or FI;


{ Writeln('        Especie       Concentracao         Numero de Mol');
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
  Writeln(especie^.esp[i].nome:15,' =',temp_nmol[i]/sistema^.vol[fase],especie^.esp[i].nmol);
  Writeln(#7);    }

  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    begin
      C:=(temp_nmol[i]/sistema^.vol[fase]);
      nmol:=temp_nmol[i];
    end;
End;

Procedure Titular( var especie :PEspecie;
                   var sistema :PSistema;
                   FI          :Boolean);
(**************************************************************************
Objetivo: Realiza a simulacao de titulacao.
Entrada :
Saida   :
Erro    :
**************************************************************************)

var                   i,n,j,np : Integer;
         Ct,Ctp,V,Vo,Vt,adnmolt,
                    adnmoltp,aux,soma : Float;
                               u : Char;
                            name : String;
                            arq  : Text;
                         sinal   : Boolean;
            { Calculo de tempo }
                       H,M,S,Sec : Array[1..2] of Word;
                       tempo     : Real;
Begin
  ClrScr;
  Writeln;
  Writeln(' Lista de especies');
  Writeln;
  For i:= 1 to especie^.nesp do
    Writeln(i,'-',especie^.esp[i].nome);
  Writeln;
  Write('Especie titulante (1 ate ',especie^.nesp,') = ');Readln(n);
  Write('Par da especie titulante (1 ate ',especie^.nesp,') =');Readln(np);
  Write('Concentracao do titulante ',especie^.esp[n].nome,' = ');Readln(Ct);
  Ctp:=Ct;
  Write('Volume do incremento na titula‡Æo (ml) = ');Readln(V);
  Write('Volume total de titulante (mL) = '); Readln(Vt);
  Writeln;
  Write('Deseja gravar os dados ? (S/N)'); u := UpCase(ReadKey);

  { 10% da aliquota do titulante }
  If (especie^.esp[n].nmol = 0.0) then especie^.esp[n].nmol:=Ct*V*1E-3*0.1;

  ClrScr;
  sinal:=true;

  { Em numero de moles }
  adnmolt := (Ct * V)*1E-3;
  adnmoltp := (Ctp * V)*1E-3;

  {Mostrar}
  Writeln('numero de mols adicionado=',adnmolt);
  Readln;

  If u = 'S' then
  begin
    sinal:=false;
    Writeln(' Gravar os dados da titulacao');
    Writeln;
    Write('Nome do arquivo:');Readln(name);
    Assign(arq,'t_'+name+'.dat');
    Rewrite(arq);

    { Gravar os dados da simulacao no arquivo }

    Writeln(arq,'Dados da simulacao');
    Writeln(arq);
    For i:= 1 to especie^.nesp do with especie^.esp[i] do
      Writeln(arq,nome,',',nmol);
    Writeln(arq,'Esp.tit.=',',',especie^.esp[n].nome);
    Writeln(arq,'Conc. do tit.=',',',Ct);
    Writeln(arq,'Vol.de incremento =',',',V,',',' mL');
    Writeln(arq,'Vol.total de titulante =',',',Vt,',',' mL');
    Writeln(arq);
  end;

  { Espalhar o sistema }
  Espalha(especie,sistema);

  ClrScr;
  aux := 0.0;
  j := 0;
  If u = 'S' then
  begin
    Write(arq,' ',',');
    For i:= 1 to especie^.nesp do with especie^.esp[i] do Write(arq,nmol,',');
    Writeln(arq);
    Write(arq,'Volume',',');
    For i:= 1 to especie^.nesp do with especie^.esp[i] do Write(arq,nome,',');
    Writeln(arq);
  end;

  { Inicio do tempo }
  GetTime(H[1],M[1],S[1],Sec[1]);

  Repeat
    Equilibrar(especie,sistema,sinal,FI);
{   Writeln;
    Write(aux:3:3,' ');
    For i:= 1 to especie^.nesp do with especie^.esp[i] do Write(nmol:14);
    Writeln;
}
    If (u = 'S') then
      begin
        Write(arq,aux:3:3,',');
        For i:= 1 to especie^.nesp do with especie^.esp[i] do Write(arq,nmol,',');
        Writeln(arq);
      end;
    Inc(j);
    aux := aux + V;
    { Efeito da diluicao }
    sistema^. vol[especie^.esp[n].fase]:=sistema^.vol[especie^.esp[n].fase]+(V*1E-3);
    Writeln('volume=', sistema^. vol[especie^.esp[n].fase]);
    especie^.esp[n].nmol := especie^.esp[n].nmol + adnmolt;
    { So' serve para efeito da forca ionica }
    {especie^.esp[np].nmol := especie^.esp[np].nmol + adnmoltp;}

  Until (aux >= Vt) ;

  { Tempo Inicial e Final na tela }
  GetTime(H[2],M[2],S[2],Sec[2]);
  tempo:= ((H[2]*3600)+(M[2]*60)+(S[2]))-((H[1]*3600)+(M[1]*60)+(S[1]));
  Writeln(H[1],':',M[1],':',S[1],':',Sec[1]);
  Writeln(H[2],':',M[2],':',S[2],':',Sec[2]);
  Write('time =',tempo:4:2);

  { Tempo Inicial e Final gravado no arquivo }
  If u = 'S' then
    begin
      Writeln(arq);
      Writeln(arq,H[1],':',M[1],':',S[1],':',Sec[1]);
      Writeln(arq,H[2],':',M[2],':',S[2],':',Sec[2]);
      Write(arq,'tempo gasto =',tempo:4:2);
      Close(arq);
    end;
End;

Procedure Distribuicao ( var especie :PEspecie;
                         var sistema :PSistema;
                             FI      :Boolean);
(**************************************************************************
Objetivo: Realiza a curva de  distribuicao das especies
Entrada :
Saida   :
Erro    :
**************************************************************************)

var
                  i,n,j,pH : Integer;
                  aux_x,aux_y : Integer;
                  aux_x1,aux_y1 : String;
         Ct,Ctp,V,Vo,Vt,adnmolt,
       incremento,aux_incremento,
         denominador,adnmoltp,aux,
  x,y,aux_pH,ph_inicial,ph_final : Float;
                         g, u, w : Char;
                            name : String;
                        arq,arq1 : Text;
                            sinal: Boolean;
            { Calculo de tempo }
                       H,M,S,Sec : Array[1..2] of Word;
                       tempo     : Real;
                       lista,lista_metal: Array[1..MAXESP] of Integer;

            GraphDriver, GraphMode: Integer;
                         code     : Integer;
Begin

 

  ClrScr;
  Writeln('Curva de distribuicao das especies em funcao do pH');
  Writeln;
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    If (nome = '(H+)') or (nome ='(H+1)') then pH:=i;
  Writeln;
  Repeat
    Writeln;
    Write('Valor do pH inicial =');Readln(ph_inicial);
    Write('Valor do pH final =');Readln(ph_final);
    Write('Incremento =');Readln(incremento);
    Writeln;
    Write('Os dados estao corretos (S/N)? ');
  Until UpCase(ReadKey) = 'S';
  Writeln;
  j:=0;
  Writeln('especie     Incluir (S/N)');
  For i:= 1 to especie^.nesp do
    begin
      Write(i,'-',especie^.esp[i].nome+'     ');u:=UpCase(ReadKey);
      If u='S' then
        begin
          Inc(j);
          lista[j]:=i;
          Write('N. de ligante / metal =');Readln(lista_metal[j]);
        end else  Writeln;
    end;
  n:=j;

  For i:= 1 to n do
    Writeln('lista[',i,'] =',lista[i],'   ',lista_metal[i]);

  Write('Deseja gravar os dados no formato de distribuicao de especies (S/N)?'); u := UpCase(ReadKey);
  If u='S' then
  begin
    Readln;
    Write('Nome  do arquivo ='); Readln(name);
    Assign(arq,{'c:\tp\equil\'+}name+'.dat');
    Rewrite(arq);
  end;
  Writeln;

  { Gravar as concentracoes das especies }
  Write('Deseja gravar os numero de mols especies (S/N)?'); w := UpCase(ReadKey);
  If w='S' then
  begin
    Readln;
    Write('Nome  do arquivo ='); Readln(name);
    Assign(arq1,name+'.dat');
    Rewrite(arq1);
  end;
  Writeln;

  Write('Deseja graficar ? (S/N)?'); g:= UpCase(ReadKey);

  { Adicionar um tampao imaginario }
  especie^.esp[especie^.nesp+1].nome:='(A)';
  especie^.esp[especie^.nesp+1].nmol:=1E-2;
  especie^.esp[especie^.nesp+1].fase:=0;
  sistema^.coef[especie^.nesp+1,sistema^.neq+1]:=1;
  especie^.esp[especie^.nesp+2].nmol:=1E-2;
  especie^.esp[especie^.nesp+2].nome:='(B)';
  especie^.esp[especie^.nesp+2].fase:=0;
  sistema^.coef[especie^.nesp+2,sistema^.neq+1]:=-1;
  especie^.nesp:=especie^.nesp+2;
  sistema^.equ[sistema^.neq+1].K:=exp(-ph_inicial*2.302585);
  sistema^.equ[sistema^.neq+1].K_1:=1/(sistema^.equ[sistema^.neq+1].K);
  sistema^.equ[sistema^.neq+1].fm:=sistema^.equ[sistema^.neq+1].K_1;
  sistema^.coef[pH,sistema^.neq+1]:=1;
  sistema^.neq:=sistema^.neq+1;

  { expo }
  sistema^.expo:=sistema^.coef;
  For i:= 1 to sistema^.neq do with sistema^ do
    For j:= 1 to especie^.nesp do
    begin
      If tf[especie^.esp[j].fase]=pura then expo[j,i]:=0;
     end;

  aux_incremento:=ph_inicial;
  sinal:=false;
  Clrscr;

  denominador:=0.0;
  For i:= 1 to n do
    denominador:=denominador+(lista_metal[i]*especie^.esp[lista[i]].nmol);

  { Inicializacao grafica }
  If (g='S') then
  begin
    DetectGraph(GraphDriver, GraphMode);
    InitGraph(GraphDriver,GraphMode,'C:\TP\BGI');
    if GraphResult <> grOk then Halt(1);
    Rectangle(10,10,629,469);
  end;

  If (u='S') then
  begin
    { Nome das especies }
    Write(arq,' '+',');
    For  i:= 1 to n do with especie^.esp[lista[i]] do
        Write(arq,nome+',');
     Writeln(arq);
  end;

  If (w='S') then
  begin
    For i:= 1 to especie^.nesp do with especie^.esp[i] do
      Writeln(arq1,nome,',',nmol,',',sistema^.vol[fase]);
    Writeln(arq1);
    For i:= 1 to (especie^.nesp-1) do with especie^.esp[i] do
      Write(arq1,nome,',');
    Writeln(arq1,especie^.esp[especie^.nesp].nome);
  end;

  { Espalhar o sistema }
  Espalha(especie,sistema);

  { Inicio da simulacao }
  GetTime(H[1],M[1],S[1],Sec[1]);

  While aux_incremento <=  ph_final do
    begin
       Equilibrar(especie,sistema,sinal,FI);

       aux_incremento:=aux_incremento+incremento;

       If sistema^.equ[sistema^.neq].fm>1 then
         sistema^.equ[sistema^.neq].fm:=(1/exp(-aux_incremento*2.302585093))
       else sistema^.equ[sistema^.neq].fm:=(exp(-aux_incremento*2.302585093));

       {Writeln('pH =',-0.434294481903*ln(especie^.esp[pH].nmol/sistema^.vol[especie^.esp[pH].fase]));}
       Writeln('pH=',aux_incremento-incremento);

       { Guardar os dados da distribuicao de especies }
       If (u = 'S') then
         begin
           Write(arq,-0.434294481903*ln(especie^.esp[pH].nmol/sistema^.vol[especie^.esp[pH].fase]):3:3,',');
           For  i:= 1 to n do with especie^.esp[lista[i]] do
             Write(arq,(lista_metal[i]*(nmol/denominador)),',');
           Writeln(arq);
         end;

       { Guardar os valores de numero de mols }
       If (w='S') then
         begin
           For i:= 1 to (especie^.nesp-1) do with especie^.esp[i] do
             Write(arq1,nmol,',');
           Writeln(arq1,especie^.esp[especie^.nesp].nmol);
         end;

        { Plotar os pontos }
        If (g='S') then
        begin
          aux_pH:=-0.434294481903*ln(especie^.esp[pH].nmol/sistema^.vol[especie^.esp[pH].fase]);
          x:=-1*(((pH_inicial-aux_pH)*(629-10))/(pH_final-pH_inicial))+10;
          Str(x:3:0,aux_x1);
          Val(aux_x1,aux_x,code);
          For i:= 1 to n do with especie^.esp[lista[i]] do
            begin
              y:=-1*(((1-(lista_metal[i]*(nmol/denominador)))*(469-10))/(0-1))+10;
              Str(y:3:0,aux_y1);
              Val(aux_y1,aux_y,code);
              PutPixel(aux_x,aux_y,3);
            end;
        end;
    end;

  If (g='S') then
  begin
    Readln;
    CloseGraph;
  end;

  { Tempo final da simulacao }
  GetTime(H[2],M[2],S[2],Sec[2]);
  tempo:= ((H[2]*3600)+(M[2]*60)+(S[2]))-((H[1]*3600)+(M[1]*60)+(S[1]));
  Writeln(H[1],':',M[1],':',S[1],':',Sec[1]);
  Writeln(H[2],':',M[2],':',S[2],':',Sec[2]);
  Write('tempo gasto =',tempo:4:2);

  { Guardar os tempos nos arquivos }
  If u = 'S' then
    begin
      Writeln(arq);
      Writeln(arq,H[1],':',M[1],':',S[1],':',Sec[1]);
      Writeln(arq,H[2],':',M[2],':',S[2],':',Sec[2]);
      Write(arq,'tempo gasto =',tempo:4:2);
      Close(arq);
    end;

  If w = 'S' then
    begin
      Writeln(arq1);
      Writeln(arq1,H[1],':',M[1],':',S[1],':',Sec[1]);
      Writeln(arq1,H[2],':',M[2],':',S[2],':',Sec[2]);
      Write(arq1,'tempo gasto =',tempo:4:2);
      Close(arq1);
    end;

End;

Procedure CalculoSimples ( var especie  :PEspecie;
                            var sistema :PSistema;
                            FI          :Boolean);
var
       w : Char;
    arq1 : Text;
    i,j  : Integer;
    name : String;
    sinal: Boolean;

    { Calculo de tempo }
    H,M,S,Sec : Array[1..2] of Word;
    tempo     : Real;
    lista,lista_metal: Array[1..MAXESP] of Integer;

Begin

  ClrScr;
  { Gravar as concentracoes das especies }
  Write('Save file numbers mols? (Y/N)?'); w := UpCase(ReadKey);
  Writeln;

  If (w='Y') then
  begin
    Write('File Name (oly name) ='); Readln(name);
    Assign(arq1,name+'.dat');
    Writeln('Save file '+name+'.dat');
    Rewrite(arq1);
  end;
  Writeln;

  If (w='Y') then
  begin
    Writeln(arq1,'Inicial concentrations and numbers mols');
    For i:= 1 to especie^.nesp do with especie^.esp[i] do
      Writeln(arq1,nome,',',nmol,',',nmol/sistema^.vol[fase]);
    Writeln(arq1);


  end;

  Clrscr;
  sinal:=false;
{ GotoXY(27,9);
  Write('Calculation ....');}
  GetTime(H[1],M[1],S[1],Sec[1]);

  Espalha(especie,sistema);

  Equilibrar(especie,sistema,sinal,FI);

  Writeln(#7);
  GetTime(H[2],M[2],S[2],Sec[2]);
  tempo:= ((H[2]*3600)+(M[2]*60)+(S[2]))-((H[1]*3600)+(M[1]*60)+(S[1]));
(*
  Writeln(H[1],':',M[1],':',S[1],':',Sec[1]);
  Writeln(H[2],':',M[2],':',S[2],':',Sec[2]);
*)
  Write('time consumed =',tempo:4:2,'s');
  GotoXY(2,24);
  Write('Next page <Enter>:');
  Readln;




  { Save file numbers mols }
  If (w='Y') then
    begin

      Writeln(arq1,'Equilibrium concentrations and numbers mols');
      For i:= 1 to especie^.nesp do with especie^.esp[i] do
        Writeln(arq1,nome,',',nmol,',',C);
      Writeln(arq1);

(*

      For i:= 1 to (especie^.nesp-1) do with especie^.esp[i] do
        Write(arq1,nmol,',');
      Writeln(arq1,especie^.esp[especie^.nesp].nmol);
*)

      Writeln(arq1);
      Writeln(arq1,H[1],':',M[1],':',S[1],':',Sec[1]);
      Writeln(arq1,H[2],':',M[2],':',S[2],':',Sec[2]);
      Write(arq1,'time comsumed =',tempo:4:2,'s');
      Close(arq1);
    end;
End;

Procedure Inicializar(var especie :PEspecie;
                      var sistema :PSistema);
(**************************************************************************
Objetivo:  Inicialzia os ponteiros
Entrada :
Saida   :
Erro    :
**************************************************************************)

Begin
  New(especie);
  New(sistema);
End;

Procedure Finalizar( var especie : PEspecie;
                     var sistema : PSistema);
(**************************************************************************
Objetivo: Finaliza os ponteiros.
Entrada :
Saida   :
Erro    :
**************************************************************************)

Begin
  Dispose(Sistema);
  Dispose(Especie);
End;

{ Programa Principal }

var
  Especie          : PEspecie;
  Sistema          : PSistema;
  nomearq,aux_nome : String;
  i,j,opcao        : Integer;
  g,r              : Char;
  sinal,FI         : Boolean;
  { Calculo de tempo }
  H,M,S,Sec : Array[1..2] of Word;
  tempo     : Real;

BEGIN
  Clrscr;
  Writeln('****************************************************************************');
  Writeln('*    SEQEx - Simulator de Equilibrios Quimicos por Extensao de reacao      *');
  Writeln('*                                                                          *');
  Writeln('*    Program Equilibrium Simulator for Multiphase Equilibri on the         *');
  Writeln('*    extent of reaction and Newton-Raphson method with    globally         *');
  Writeln('*    convegent strategy                                                    *');
  Writeln('*                                                                          *');
  Writeln('*    Universidade de Sao Paulo - Instituto de Quimica - Sao Paulo - Brazil *');
  Writeln('*                                                                          *');
  Writeln('*    Created by  Robinson Mioshi - rnmioshi@quim.iq.usp.br                 *');
  Writeln('*                Claudimir L. do Lago - claudemi@quim.iq.usp.br            *');
  Writeln('*                                                                          *');
  Writeln('*    Version 1.1.0                                                         *');
  Writeln('*    Created on 30/11/94                   Last Modified: 17/02/97         *');
  Writeln('*                                                                          *');
  Writeln('*    Language Program: Borland Pascal 6.0 (DOS)                            *');
  Writeln('****************************************************************************');

  GotoXY(2,24);
  Write(' Next Page <Enter> :');
  Readln;

  Clrscr;



  Inicializar(especie,sistema);




  Write('Name data file (only name) =');
  Readln(aux_nome);

  { Data input file data equilibrium expressions}
  (* Example:  {H2O} = {H+} +  {OH-}            Kw = 1.8E-16
               {CO2}^1 + {H2O} = {CO2aq}        K = 2.55E-5
               {CO2aq} =  {H+} + {HCO3-}        K = 4.45E-7
               {HCO3-} = {H+} + {CO3-2}         K = 4.69E-11
               {Ca+2} + {CO3-2} = {CaCO3}^2     K = 2.239E8   *)
  nomearq:='d_'+aux_nome+'.dat';
  LeReacao(nomearq,sistema,especie);
  Writeln('Read d_'+aux_nome+'.dat');



  { Data input file numbers mols and effective diameter of the hydrated }
  (* Example:  {H2O} = 55.6E-2     0.0
             {CO2}^1 = 7.95E6      0.0
             {H+}    = 1E-9        0.0
             {OH-}   = 1E-9        3.5
             {CO2aq} = 1.0E-100    0.0
             {Ca+2}  = 1.0E-70     6.0
             {CaCO3}^2 = 1.0       0.0                      *)
  nomearq:='c_'+aux_nome+'.dat';
  LeConcentr(nomearq,especie);
  Writeln('Read c_'+aux_nome+'.dat');

  { Data input file volume multicomponent phase and volume}
  (* Example:  M0 = 10E-3
               M1 = 1E9
               P2                                           *)
  nomearq:='v_'+aux_nome+'.dat';
  LeVolume(nomearq,sistema,especie);
  Writeln('Read v_'+aux_nome+'.dat');
  Writeln;

  Write('Include ionic strength effect? (Y/N)'); r := UpCase(ReadKey);
  Writeln;
  If (r='Y') then
  begin
    FI:=false;

    (*
    Writeln('Effective diameter of the hydrated (angstrom)');
    For i:= 1 to especie^.nesp do with especie^.esp[i] do
      begin
        If carga<>0 then
          begin
            Writeln(nome,'  =',diametro:2);
          end;
      end;
    GotoXY(2,24);
    Write('Next page <Enter> :');
    Readln;       *)
  end else FI:=true;

  { Imprime  nomes, cargas e fases das especies }
  Writeln;
  Clrscr;
  Writeln('Components':15,'    Numbers mols':6,'  Volume':6,'        Diameter':6);
  For i:= 1 to especie^.nesp do with especie^.esp[i] do
    Writeln(nome:15,'   ',nmol:6,'    ',sistema^.vol[fase]:6,'    ',diametro:6);
  GotoXY(2,24);
  Write('Next page <Enter>:');
  Readln;

(*
  { Matriz dos coeficientes }
  Writeln;
  Clrscr;
  Writeln('Coluna : Especies');
  Writeln('Linha  : Reacoes');
  For i:= 1 to sistema^.neq do with sistema^ do
    begin
      For j:= 1 to especie^.nesp do with especie^ do
        Write(coef[j,i]:2,' ');
      Writeln;
    end;
  GotoXY(2,24);
  Write(' Aperte qualquer tecla pra continuar <Enter> :');
  Readln;

  Clrscr;

  Writeln('Espalhar as concentra‡äes iniciais das esp‚cies');
  Writeln;
  Writeln('        Esp‚cie         Concentra‡Æo           N£mero de Mol');
  For i:= 1 to especie^.nesp do
    Writeln(especie^.esp[i].nome:15,' ',(especie^.esp[i].nmol/sistema^.vol[especie^.esp[i].fase])
            ,especie^.esp[i].nmol);
  Writeln;
  Write(' Aperte qualquer tecla pra continuar <Enter> :');
  Readln;

  Clrscr;
  Writeln('Escolha a opcao:');
  Writeln;
  Writeln('1 - Calcular as concentracoes');
  Writeln('2 - Titulacao');
  Writeln('3 - Curva de distribuicao das especies');
  Writeln;
  Write('Opcao (1,2 ou 3) :');Readln(opcao);

  Case opcao of
    1: CalculoSimples(especie,sistema,FI);
    2: Titular(especie,sistema,FI);
    3: Distribuicao(especie,sistema,FI);
  end;

*)
  CalculoSimples(especie,sistema,FI);

    { ******************* Saida de Tela ********************* }

    Clrscr;
    Writeln;
    Writeln('Equilibrium concentration components'); Writeln;
    Writeln('Components':15,'  Concentration (mol/L)':15,'Numbers mols':15);
    For i:= 1 to especie^.nesp do
      Writeln(especie^.esp[i].nome:15,' =',especie^.esp[i].C,especie^.esp[i].nmol);
    Readln;

   Finalizar(especie,sistema);

END.

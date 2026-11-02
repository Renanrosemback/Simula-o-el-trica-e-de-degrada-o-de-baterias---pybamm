# Otimização de Baterias - BESS

## Sobre o projeto

Este projeto consiste no desenvolvimento de um modelo computacional para a simulação e otimização de um sistema de armazenamento de energia por baterias (BESS — Battery Energy Storage System), aplicado ao sistema de tração de um caminhão fora de estrada utilizado em operações de mineração.

O modelo realiza a simulação do comportamento elétrico da bateria ao longo de um horizonte de análise, considerando o perfil de potência do sistema, o estado de carga (SOC), os limites de operação e os efeitos de degradação da bateria. A partir dos resultados da simulação, também são avaliados aspectos energéticos e econômicos do sistema, incluindo o custo de ciclo de vida (LCC — Life Cycle Cost).

O projeto foi desenvolvido em Python e estruturado utilizando programação orientada a objetos (POO), buscando separar as diferentes responsabilidades do modelo em classes e módulos independentes.

## Objetivos

Os principais objetivos do projeto são:

- Desenvolver um modelo computacional em Python capaz de representar o comportamento elétrico de uma bateria submetida a um perfil de potência de operação;

- Simular a evolução do estado de carga (SOC) e dos parâmetros elétricos da bateria ao longo do tempo;

- Considerar os efeitos da degradação cíclica e da degradação de calendário sobre o desempenho da bateria;

- Avaliar indicadores energéticos e econômicos do sistema, incluindo a energia processada, a energia regenerativa aproveitada e o custo de ciclo de vida (LCC);

- Determinar, por meio de um processo de otimização, uma configuração adequada para o número de strings de baterias conectadas em paralelo (Np), buscando minimizar o custo de ciclo de vida do sistema.

## Estrutura do projeto

O projeto foi desenvolvido de forma modular, utilizando programação orientada a objetos (POO) para organizar as diferentes etapas envolvidas na modelagem, simulação, análise econômica e otimização do sistema de armazenamento de energia. A estrutura foi pensada de forma que cada módulo possua uma responsabilidade específica, permitindo que os diferentes componentes do modelo sejam desenvolvidos e modificados de maneira independente, enquanto as classes responsáveis pela simulação fazem a integração entre eles.

De maneira geral, o funcionamento do projeto pode ser entendido a partir de uma sequência de etapas. Inicialmente, o sistema recebe um perfil de potência que representa a demanda do caminhão ao longo de sua operação. Esse perfil é tratado e convertido para o formato utilizado pelo modelo da bateria. A partir dele, a classe responsável pela simulação coordena o comportamento elétrico da bateria, calcula a evolução do estado de carga (SOC) e determina a potência efetivamente processada pelo sistema de armazenamento. Ao final de cada período de simulação, os resultados obtidos são utilizados pelo modelo de degradação para atualizar o estado da bateria. Esse processo é repetido ao longo do horizonte de análise.

A organização do projeto pode ser dividida em cinco componentes principais: configuração dos parâmetros, modelagem elétrica da bateria, modelagem da degradação, simulação do sistema e análise econômica e otimização.

### Configuração dos parâmetros

O arquivo `config.py` concentra os principais parâmetros e constantes utilizados pelo modelo. A centralização desses valores evita que os mesmos parâmetros precisem ser definidos repetidamente em diferentes arquivos e facilita a alteração das condições de simulação.

Entre os parâmetros utilizados estão as características elétricas da bateria, limites de estado de carga, capacidade nominal, taxa máxima de corrente, número de células em série, número de strings em paralelo, parâmetros relacionados aos custos do sistema, temperatura de operação e horizonte de análise.

Dessa forma, os demais módulos podem utilizar uma configuração comum, mantendo uma maior consistência entre as diferentes etapas do projeto.

### Modelagem elétrica da bateria

A modelagem elétrica é implementada no arquivo `battery_circuit_POO.py`, que contém a classe `Bateria`. Essa classe representa o comportamento elétrico do banco de baterias durante a simulação.

O modelo utiliza um circuito elétrico equivalente composto por uma tensão de circuito aberto dependente do estado de carga, uma resistência série e dois ramos ZARC. Cada ramo ZARC é representado por cinco sub-ramos RC, permitindo representar diferentes características dinâmicas da resposta da bateria.

A classe `Bateria` é responsável por receber a potência solicitada pelo sistema e determinar a corrente correspondente, considerando a tensão instantânea da bateria e os limites de corrente definidos para sua operação. A partir da corrente calculada, o modelo determina a tensão do banco de baterias e atualiza o estado de carga por meio de contagem de Coulomb.

Durante a simulação, também são considerados os limites mínimo e máximo de SOC. Dessa forma, a potência efetivamente processada pela bateria pode ser limitada quando o estado de carga atinge uma condição que impede a continuidade da operação.

A classe também possui métodos relacionados à atualização dos parâmetros elétricos da bateria. Esses parâmetros podem ser alterados em função da degradação acumulada ao longo da vida útil do sistema. Assim, a modelagem elétrica não representa apenas uma bateria ideal ou estática, mas uma bateria cujas características elétricas podem evoluir ao longo do horizonte de análise.

### Modelagem da degradação

O envelhecimento da bateria é tratado separadamente da modelagem elétrica, por meio da classe `BatteryDegradation`, implementada no arquivo `lifetime_degradationPOO.py`.

A degradação é dividida em dois mecanismos principais: degradação cíclica e degradação de calendário. A degradação cíclica está relacionada ao uso da bateria e ao comportamento do SOC durante os ciclos de operação, enquanto a degradação de calendário representa o envelhecimento associado à passagem do tempo e às condições de operação.

A classe recebe os resultados da simulação do SOC e utiliza essas informações para determinar a evolução da degradação ao longo dos meses. Os resultados são utilizados para calcular a perda de capacidade da bateria e as alterações de seus parâmetros elétricos.

Essa separação entre o modelo elétrico e o modelo de degradação permite que cada componente tenha uma responsabilidade bem definida. A classe `Bateria` representa o comportamento elétrico instantâneo do sistema, enquanto `BatteryDegradation` determina como as características da bateria se modificam ao longo de sua vida útil.

### Perfil de potência

O tratamento do perfil de potência é realizado pela classe `PerfilPotencia`, presente no arquivo `SOCestimation_testePOO.py`.

Essa classe é responsável por carregar os dados de potência utilizados como entrada da simulação e realizar o tratamento necessário para adequá-los ao modelo. O perfil contém informações relacionadas à potência de tração e à potência de frenagem regenerativa do caminhão.

A partir desses dados, é obtido o perfil de potência que será utilizado pela bateria durante a simulação. O perfil é expandido para representar o período de análise definido para cada mês e convertido para as unidades e convenções de sinais utilizadas pelo modelo elétrico.

A existência de uma classe específica para o perfil de potência permite separar o tratamento dos dados de entrada da lógica responsável pela simulação da bateria. Dessa maneira, diferentes perfis de operação podem ser utilizados sem que seja necessário modificar diretamente o modelo elétrico.

### Simulação

A integração entre os diferentes componentes ocorre principalmente por meio da classe `Simulacao`, também implementada em `SOCestimation_testePOO.py`.

A classe `Simulacao` funciona como o elemento responsável por coordenar a execução do modelo ao longo do horizonte de análise. Para isso, ela possui objetos das classes `Bateria` e `BatteryDegradation`, além de utilizar um objeto `PerfilPotencia` para fornecer o perfil de operação do sistema.

Durante cada mês de simulação, a classe executa inicialmente o modelo elétrico da bateria utilizando o perfil de potência. Como resultado, são obtidos o estado de carga ao longo do período e a potência efetivamente processada pela bateria.

Em seguida, os resultados energéticos são calculados, incluindo a energia processada pelo banco de baterias e a energia regenerativa efetivamente aproveitada. O histórico de SOC obtido durante a simulação é então enviado ao modelo de degradação.

Após o cálculo da degradação, os parâmetros da bateria são atualizados para representar o estado degradado ao final daquele período. Esse processo permite que a simulação do mês seguinte seja realizada utilizando as características atualizadas da bateria.

A classe `Simulacao` também controla a ocorrência de substituições do banco de baterias. Quando a degradação acumulada ultrapassa o limite de fim de vida estabelecido, o banco pode ser substituído e seus parâmetros são reinicializados. Dessa forma, a simulação consegue representar não apenas o envelhecimento contínuo da bateria, mas também as possíveis substituições ao longo do horizonte de análise.

O fluxo geral da simulação pode ser representado de forma simplificada como:

`Perfil de potência → Bateria → SOC e potência → Degradação → Atualização da bateria → Próximo mês`

Esse ciclo é repetido até que todos os meses definidos no horizonte de análise tenham sido simulados.

### Análise do custo de ciclo de vida

Após a estrutura de simulação estar estabelecida, os resultados podem ser utilizados para realizar a análise econômica do sistema. Essa etapa é implementada no arquivo `lcc.py`.

O módulo de LCC (Life Cycle Cost) calcula o custo total associado à utilização do sistema de armazenamento ao longo do horizonte de análise. Para isso, são considerados diferentes componentes econômicos, como o investimento inicial no banco de baterias e no conversor, os custos associados ao consumo de energia, os benefícios relacionados à energia regenerativa aproveitada, os custos de substituição da bateria e o valor residual do banco ao final do período analisado.

Um dos principais parâmetros avaliados pelo modelo econômico é o número de strings conectadas em paralelo (`Np`). Esse parâmetro influencia diretamente a capacidade instalada do banco de baterias e, consequentemente, tanto o investimento inicial quanto o comportamento operacional e a degradação do sistema.

Para evitar que o perfil de potência precise ser carregado novamente a cada avaliação, o módulo também possui uma estrutura de armazenamento temporário dos dados do perfil, permitindo que diferentes configurações de `Np` sejam avaliadas utilizando os mesmos dados de entrada.

Assim, o módulo de LCC funciona como uma camada entre a simulação física do sistema e o processo de otimização. Para cada configuração avaliada, a simulação fornece os resultados operacionais e de degradação, que são posteriormente convertidos em indicadores econômicos.

### Otimização

O processo de otimização é implementado no arquivo `otimizacao.py`. Seu objetivo é determinar o número de strings em paralelo que proporciona o menor custo de ciclo de vida para o sistema.

Para isso, o processo de otimização utiliza uma estratégia em duas etapas. Inicialmente, é realizada uma busca em uma grade de diferentes valores de `Np`, permitindo identificar uma região próxima à configuração de menor custo. Em seguida, essa região é refinada utilizando um método de otimização escalar baseado em Brent.

Para cada valor de `Np` avaliado, o módulo de otimização solicita ao modelo de LCC o cálculo do custo de ciclo de vida. O LCC, por sua vez, executa a simulação do sistema para aquela configuração e utiliza os resultados obtidos para calcular o custo total.

Dessa forma, existe uma relação hierárquica entre os diferentes componentes do projeto:

`Otimização → LCC → Simulação → Bateria + Degradação + Perfil de potência`

A otimização representa, portanto, a camada mais externa do modelo. Ela não realiza diretamente os cálculos elétricos ou de degradação, mas utiliza os resultados produzidos pelas camadas inferiores para avaliar diferentes configurações do sistema.

### Integração entre os módulos

A principal característica da estrutura do projeto é a separação das responsabilidades entre os diferentes componentes. O modelo elétrico da bateria não precisa conhecer os detalhes do cálculo econômico, assim como o módulo de otimização não precisa conhecer diretamente as equações utilizadas para calcular a tensão ou o SOC da bateria.

Essa organização permite que o projeto seja dividido em diferentes níveis de abstração. Em um nível mais baixo está o modelo da bateria, responsável pela representação do comportamento elétrico. Em seguida está o modelo de degradação, responsável pela evolução das características da bateria ao longo do tempo. A classe `Simulacao` integra esses modelos e reproduz a operação do sistema durante o horizonte de análise. Acima da simulação encontra-se o cálculo do custo de ciclo de vida, que transforma os resultados físicos em indicadores econômicos. Finalmente, o módulo de otimização utiliza o LCC para comparar diferentes configurações e determinar a configuração mais adequada segundo o critério de custo adotado.

Essa estrutura também facilita a manutenção e a expansão do projeto. Alterações no modelo elétrico podem ser realizadas na classe `Bateria` sem que seja necessário modificar diretamente o algoritmo de otimização. Da mesma forma, mudanças no modelo de degradação podem ser incorporadas à classe `BatteryDegradation` sem alterar a lógica principal da análise econômica.

Portanto, a utilização da programação orientada a objetos permite representar o sistema de maneira modular, em que cada classe possui uma função específica e as diferentes partes do projeto se comunicam por meio de suas interfaces. O resultado é uma estrutura organizada em camadas, na qual os modelos físicos, a simulação, a análise econômica e a otimização são integrados para representar o problema completo de dimensionamento do BESS.

## Modelo do sistema

## Funcionamento

## Requisitos

## Instalação

## Como executar

## Parâmetros de configuração

## Resultados

## Referências

## Autores
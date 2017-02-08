>> daq.getDevices()

ans = 

ni: National Instruments PCIe-6351 (Device ID: 'Dev1')
   Analog input subsystem supports:
      7 ranges supported
      Rates from 0.1 to 1250000.0 scans/sec
      16 channels ('ai0' - 'ai15')
      'Voltage' measurement type
   
   Analog output subsystem supports:
      -5.0 to +5.0 Volts,-10 to +10 Volts ranges
      Rates from 0.1 to 2857142.9 scans/sec
      2 channels ('ao0','ao1')
      'Voltage' measurement type
   
   Digital subsystem supports:
      Rates from 0.1 to 10000000.0 scans/sec
      24 channels ('port0/line0' - 'port2/line7')
      'InputOnly','OutputOnly','Bidirectional' measurement types
   
   Counter input subsystem supports:
      Rates from 0.1 to 100000000.0 scans/sec
      4 channels ('ctr0','ctr1','ctr2','ctr3')
      'EdgeCount','PulseWidth','Frequency','Position' measurement types
   
   Counter output subsystem supports:
      Rates from 0.1 to 100000000.0 scans/sec
      4 channels ('ctr0','ctr1','ctr2','ctr3')
      'PulseGeneration' measurement type
   
                       AutoSyncDSA: false
                       NumberOfScans: 500000
                   DurationInSeconds: 1
                                Rate: 500000
                        IsContinuous: false
      NotifyWhenDataAvailableExceeds: 250000
IsNotifyWhenDataAvailableExceedsAuto: false
          NotifyWhenScansQueuedBelow: 250000
    IsNotifyWhenScansQueuedBelowAuto: true
              ExternalTriggerTimeout: 10
                      TriggersPerRun: 1
                              Vendor: National Instruments
                            Channels: [1x2 daq.ni.AnalogInputVoltageChannel]
                         Connections: [1x1 daq.ni.StartTriggerConnection]
                           IsRunning: false
                           IsLogging: false
                              IsDone: false
         IsWaitingForExternalTrigger: false
                   TriggersRemaining: 1
                           RateLimit: [0.1 500000.0]
                         ScansQueued: 0
               ScansOutputByHardware: 0
                       ScansAcquired: 0
							  
Analog input subsystem supports:
   7 ranges supported
   Rates from 0.1 to 1250000.0 scans/sec
   16 channels ('ai0' - 'ai15')
   'Voltage' measurement type
Properties, Methods, Events

Analog output subsystem supports:
   -5.0 to +5.0 Volts,-10 to +10 Volts ranges
   Rates from 0.1 to 2857142.9 scans/sec
   2 channels ('ao0','ao1')
   'Voltage' measurement type
Properties, Methods, Events

Digital subsystem supports:
   Rates from 0.1 to 10000000.0 scans/sec
   24 channels ('port0/line0' - 'port2/line7')
   'InputOnly','OutputOnly','Bidirectional' measurement types
Properties, Methods, Events

Counter input subsystem supports:
   Rates from 0.1 to 100000000.0 scans/sec
   4 channels ('ctr0','ctr1','ctr2','ctr3')
   'EdgeCount','PulseWidth','Frequency','Position' measurement types
Properties, Methods, Events

Counter output subsystem supports:
   Rates from 0.1 to 100000000.0 scans/sec
   4 channels ('ctr0','ctr1','ctr2','ctr3')
   'PulseGeneration' measurement type
	
     IsSimulated: false
       Terminals: [73x1 cell]
          Vendor: National Instruments
              ID: 'Dev1'
           Model: 'PCIe-6351'
      Subsystems: [1x5 daq.ClockedSubsystemInfo]
     Description: 'National Instruments PCIe-6351'
RecognizedDevice: true
d.Terminals

ans =

  73×1 cell array

    'Dev1/PFI0' %#ok<*NOPTS>
    'Dev1/PFI1'
    'Dev1/PFI2'
    'Dev1/PFI3'
    'Dev1/PFI4'
    'Dev1/PFI5'
    'Dev1/PFI6'
    'Dev1/PFI7'
    'Dev1/PFI8'
    'Dev1/PFI9'
    'Dev1/PFI10'
    'Dev1/PFI11'
    'Dev1/PFI12'
    'Dev1/PFI13'
    'Dev1/PFI14'
    'Dev1/PFI15'
    'Dev1/APFI0'
    'Dev1/RTSI0'
    'Dev1/RTSI1'
    'Dev1/RTSI2'
    'Dev1/RTSI3'
    'Dev1/RTSI4'
    'Dev1/RTSI5'
    'Dev1/RTSI6'
    'Dev1/RTSI7'
    'Dev1/20MHzTimebase'
    'Dev1/100MHzTimebase'
    'Dev1/10MHzRefClock'
    'Dev1/ChangeDetectionEvent'
    'Dev1/WatchdogExpiredEvent'
    'Dev1/WatchdogExpirationTrigger'
    'Dev1/AnalogComparisonEvent'
    'Dev1/100kHzTimebase'
    'Dev1/None'
    'Dev1/Ctr0Source'
    'Dev1/Ctr1Source'
    'Dev1/Ctr2Source'
    'Dev1/Ctr3Source'
    'Dev1/Ctr0Gate'
    'Dev1/Ctr1Gate'
    'Dev1/Ctr2Gate'
    'Dev1/Ctr3Gate'
    'Dev1/Ctr0Aux'
    'Dev1/Ctr1Aux'
    'Dev1/Ctr2Aux'
    'Dev1/Ctr3Aux'
    'Dev1/Ctr0SampleClock'
    'Dev1/Ctr1SampleClock'
    'Dev1/Ctr2SampleClock'
    'Dev1/Ctr3SampleClock'
    'Dev1/Ctr0ArmStartTrigger'
    'Dev1/Ctr1ArmStartTrigger'
    'Dev1/Ctr2ArmStartTrigger'
    'Dev1/Ctr3ArmStartTrigger'
    'Dev1/Ctr0InternalOutput'
    'Dev1/Ctr1InternalOutput'
    'Dev1/Ctr2InternalOutput'
    'Dev1/Ctr3InternalOutput'
    'Dev1/Ctr0A'
    'Dev1/Ctr1A'
    'Dev1/Ctr2A'
    'Dev1/Ctr3A'
    'Dev1/Ctr0B'
    'Dev1/Ctr1B'
    'Dev1/Ctr2B'
    'Dev1/Ctr3B'
    'Dev1/Ctr0Z'
    'Dev1/Ctr1Z'
    'Dev1/Ctr2Z'
    'Dev1/Ctr3Z'
    'Dev1/PairedCtrInternalOutput'
    'Dev1/PairedCtrOutputPulse'
    'Dev1/FrequencyOutput'

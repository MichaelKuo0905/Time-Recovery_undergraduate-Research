//==============================================================================
// ASIX AX88179/178A Windows 7 driver
//
// This document describes the driver's configurable parameters.
//==============================================================================

1. SpeedDuplex: Set the Ethernet link speed.
   -Auto:         Auto detect the Ethernet link speed
   -10Half:       Set the Ethernet force on 10HD
   -10Full:       Set the Ethernet force on 10FD
   -100Half:      Set the Ethernet force on 100HD
   -100Full:      Set the Ethernet force on 100FD
   -1000Full:     Set the Ethernet force on 1000FD


2. NetworkAddress : Allows user to set a MAC address of the device or use the 
                    device default address.

3. FlowControl: Configure flow control advertised capabilities
   -Disabled       Disable flow control
   -TxEnabled      Enable transmit flow control
   -RxEnabled      Enable receive flow control
   -RxTxEnabled    Enable transmit and receive flow control

4. WakeOnLinkChange: Wake up the computer when device detects Ehernet link Changed
   -Disabled       Disable WakeOnLinkChange function
   -Enabled        Enable WakeOnLinkChange function

5. Packet Priority & VLAN: Enable or disable the ability to insert the 802.1Q
                           priority and VLAN tags into the transmit packets.
   -Packet Priority & VLAN Disabled  Disable to insert the 802.1Q priority and VLAN tag
   -Packet Priority Enabled          Only enable to insert the 802.1Q priority tag
   -VLAN Enabled                     Only enable to insert the 802.1Q VLAN tag
   -Packet Priority & VLAN Enabled   Enable to insert the 802.1Q priority and VLAN tag

6. VLAN ID: If user set a valid VLAN ID, the driver inserts the VLAN tag with 
            this VLAN ID into the transmit packets and device will filter the
            received packets.

7. WakeOnMagicPacket: Wake up the computer when device receives a Magic Packet
   -Disabled       Disable WakeOnMagicPacket function
   -Enabled        Enable WakeOnMagicPacket function

8. WakeOnPattern: Wake up the computer when device receives a packet that matches a
                  specified pattern
   -Disabled       Disable WakeOnPattern function
   -Enabled        Enable WakeOnPattern function

9. TCPChecksumOffloadV4: Enable or disable the device to calculate the
                         TCP checksum of the transmit IPv4 packets and
                         check the TCP checksum of the received IPv4
                         packets.
   -Disabled       Disable the TCP Checksum Offload
   -TxEnabled      Enable the TCP Checksum Offload for transmit packets
   -RxEnabled      Enable the TCP Checksum Offload for received packets
   -RxTxEnabled    Enable the TCP Checksum Offload for transmit and received packets

10. UDPChecksumOffloadV4: Enable or disable the device to calculate the
                         UDP checksum of the transmit IPv4 packets and
                         check the UDP checksum of the received IPv4 packets.
   -Disabled       Disable the UDP Checksum Offload
   -TxEnabled      Enable the UDP Checksum Offload for transmit packets
   -RxEnabled      Enable the UDP Checksum Offload for received packets
   -RxTxEnabled    Enable the UDP Checksum Offload for transmit and received packets

11. TCPChecksumOffloadV6: Enable or disable the device to calculate the
                          TCP checksum of the transmit IPv6 packets and
                          check the TCP checksum of the received IPv6 packets.
   -Disabled       Disable the TCP Checksum Offload
   -TxEnabled      Enable the TCP Checksum Offload for transmit packets
   -RxEnabled      Enable the TCP Checksum Offload for received packets
   -RxTxEnabled    Enable the TCP Checksum Offload for transmit and received packets

12. UDPChecksumOffloadV6: Enable or disable the device to calculate the
                          UDP checksum of the transmit IPv6 packets and
                          check the UDP checksum of the received IPv6 packets.
   -Disabled       Disable the UDP Checksum Offload
   -TxEnabled      Enable the UDP Checksum Offload for transmit packets
   -RxEnabled      Enable the UDP Checksum Offload for received packets
   -RxTxEnabled    Enable the UDP Checksum Offload for transmit and received packets

13. IPChecksumOffloadV4: Enable or disable the device to calculate the IP 
                         checksum of the transmit IPv4 packets and check the
                         IP checksum of the received IPv4 packets.
    -Disabled      Disable the IP Checksum Offload
    -TxEnabled     Enable the IP Checksum Offload for transmit packets
    -RxEnabled     Enable the IP Checksum Offload for received packets
    -RxTxEnabled   Enable the IP Checksum Offload for transmit and received packets

14. AutoDetach: Enable or disable AutoDetach ability. if AutoDetach is enabled, 
                3 seconds later after Ethernet cable was unpluged, the device will
                be detached from USB.
    -Disabled      Disable AutoDetach
    -Enabled       Enable AutoDetach
    -Use EEPROM    Setting Disable or enable AutoDetach accoring to the EEPROM setting

15. MaskTimer: If wake up ability is enabled, the wake up function will delay this
               time to active.
    -0,4,8,12,16,20,24,28 seconds

16. WolLinkSpeed: Set the Ethernet link speed when device sleeps if wake up 
                  ability is enabled.
    -Enabled       The Ethernet link speed works on 10 Mbps first if available
    -Disabled      The Ethernet link speed works on maximum link speed

17. SelectiveSuspend: Allows NDIS to suspend an idle AX88772B network adapter by
                     transitioning the adapter to a low-power state
   -Disabled       Disable SelectiveSuspend function
   -Enabled        Enable SelectiveSuspend function

18. SSIdleTimeout: Selective suspend idle time-out in units of seconds

19. ArpOffload: When enable this ability, the device will reply the ARP request
                packet when computer is suspending. This ability is activated only
                if WOL is enabled.
    -Disabled      Disable ARP Offload
    -Enabled       Enable ARP Offload

20. NSOffload: When enable this ability, the device will reply the neighbor
               solicitation packet when computer is suspending. This ability is 
               activated only if WOL is enabled.
    -Disabled      Disable NS Offload
    -Enabled       Enable NS Offload

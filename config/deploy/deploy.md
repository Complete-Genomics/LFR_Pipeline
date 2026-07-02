
sudo pacman -S cifs-utils smbclient

smbclient -L //SERVERNAME -U DOMAIN\\username

sudo mkdir -p /mnt/storage

sudo mount -t cifs //SERVERNAME/SHARENAME /mnt/storage \
-o username=username,domain=DOMAIN,vers=3.0

sudo mount -t cifs //SERVERNAME/SHARENAME /mnt/storage \
-o credentials=/root/.smbcred,vers=3.0


smbclient -L //172.23.143.7 -U ycai

ip addr | grep inet
    inet 127.0.0.1/8 scope host lo
    inet6 ::1/128 scope host noprefixroute
    inet 172.23.143.7/24 brd 172.23.143.255 scope global dynamic noprefixroute eno1
    inet6 fe80::bfa9:ae6c:369f:401d/64 scope link
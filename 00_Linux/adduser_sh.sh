sudo adduser username -g gbh0319
sudo cp /home/ec2-user/.bashrc /home/username/.bashrc
sudo cp /home/ec2-user/.bash_profile /home/username/.bash_profile
sudo su - username -c "mkdir /home/username/.ssh"
sudo su - username -c "chmod 700 /home/username/.ssh
"
sudo su - username -c "ssh-keygen -y -f /home/ec2-user/info/username.pem > /home/username/.ssh/authorized_keys"
sudo su - username -c "ssh-keygen -y -f /home/ec2-user/info/username.pem"
ls
chmod 444 username.pem 
sudo su - username -c "ssh-keygen -y -f /home/ec2-user/info/username.pem"
sudo su - username -c "ssh-keygen -y -f /home/ec2-user/info/username.pem > /home/username/.ssh/authorized_keys"
sudo su - username -c "chmod 600 /home/username/.ssh/authorized_keys"
sudo su - username -c "source /home/username/.bashrc; source /home/username/.bash_profile"
sudo su - username -c "conda init"
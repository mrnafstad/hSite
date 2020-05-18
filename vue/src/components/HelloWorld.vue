<template>
  <div>
    <img id="profilepic" src="../assets/profilepic.jpg"> <br/>
    <h>{{ greeting }}!</h>
    <AboutMe />
  </div>
  
</template>

<script>
import AboutMe from './AboutMe.vue'
export default {
  components: {
    AboutMe
  },
    data: function() {
      return {
        greeting: "Hello world, it's me"
      };
    }
  };
</script>

<style scoped>
  h {
    padding-top: 10px;
    font-size: 2em;
    text-align: center;
    color: red;
  }
  #profilepic {
    border-radius: 50%;
    width: 300px;
    padding: 10px;
}

</style>